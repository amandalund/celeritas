//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file PhotoelectricInteractor.i.hh
//---------------------------------------------------------------------------//

#include "base/ArrayUtils.hh"
#include "random/distributions/UniformRealDistribution.hh"
#include "MockXsCalculator.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with shared and state data.
 */
CELER_FUNCTION PhotoelectricInteractor::PhotoelectricInteractor(
    const PhotoelectricInteractorPointers& shared,
    const LivermoreParamsPointers&         data,
    const ParticleTrackView&               particle,
    const Real3&                           inc_direction,
    SecondaryAllocatorView&                allocate)
    : shared_(shared)
    , data_(data)
    , inc_direction_(inc_direction)
    , inc_energy_(particle.energy().value())
    , allocate_(allocate)
    , calc_micro_xs_(shared, data, particle)
{
    REQUIRE(inc_energy_ > this->min_incident_energy()
            && inc_energy_ <= this->max_incident_energy());
    REQUIRE(particle.def_id() == shared_.gamma_id);

    inv_energy_ = 1. / inc_energy_.value();
}

//---------------------------------------------------------------------------//
/*!
 * Sample using the Livermore model for the photoelectric effect.
 */
template<class Engine>
CELER_FUNCTION Interaction
PhotoelectricInteractor::operator()(Engine& rng, ElementDefId el_id)
{
    // Select target atom
    // ElementSelector    select_element(mat_, calc_micro_xs_);
    // ElementComponentId component_id = select_element(rng);
    // ElementDefId el_id
    //    = mat_.material_view().elements()[component_id.get()].element;

    // Allocate space for the single electron to be emitted
    Secondary* photoelectron = this->allocate_(1);
    if (photoelectron == nullptr)
    {
        // Failed to allocate space for a secondary
        return Interaction::from_failure();
    }

    // Get the cross section data for the sampled element
    const LivermoreElement& el = data_.elements[el_id.get()];

    // Sample the shell from which the photoelectron is emitted
    // TODO: don't use energy = max(energy, min binding energy) here?
    real_type cutoff   = generate_canonical(rng) * calc_micro_xs_(el_id);
    real_type xs       = 0.;
    size_type shell_id = 0;
    do
    {
        const auto& shell = el.shells[shell_id++];
        if (inc_energy_ > shell.binding_energy)
        {
            if (inc_energy_ < el.thresh_low)
            {
                // Use the tabulated subshell cross sections
                XsCalculator calc_xs(shell.xs);
                xs += inv_energy_ * calc_xs(inc_energy_.value());
            }
            else
            {
                // Use parameterized integrated subshell cross sections
                const auto& param = inc_energy_ >= el.thresh_high
                                        ? shell.param_high
                                        : shell.param_low;

                // Calculate the subshell cross section from the fit parameters
                // and energy as \sigma(E) = a_1 / E + a_2 / E^2 + a_3 / E^3 +
                // a_4 / E^4 + a_5 / E^5 + a_6 / E^6.
                // clang-format off
                xs = inv_energy_ * (param[0] + inv_energy_ * (param[1]
                   + inv_energy_ * (param[2] + inv_energy_ * (param[3]
                   + inv_energy_ * (param[4] + inv_energy_ * param[5])))));
                // clang-format on
            }
        }
    } while (xs < cutoff && shell_id != el.shells.size() - 1);

    // Construct interaction for change to primary (incident) particle
    Interaction result = Interaction::from_absorption();

    // If the binding energy of the sampled shell is greater than the incident
    // photon energy, don't produce secondaries
    if (el.shells[shell_id].binding_energy > inc_energy_)
    {
        return result;
    }

    // Outgoing secondary is an electron
    result.secondaries    = {photoelectron, 1};
    photoelectron->def_id = shared_.electron_id;

    // Electron kinetic energy is the difference between the incident photon
    // energy and the binding energy of the shell
    photoelectron->energy = MevEnergy{
        inc_energy_.value() - el.shells[shell_id].binding_energy.value()};

    // Direction of the emitted photoelectron is sampled from the
    // Sauter-Gavrila distribution
    photoelectron->direction = this->sample_direction(rng);

    // TODO: Atomic relaxation

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Sample a direction according to the Sauter-Gavrila distribution.
 *
 * \note The Sauter-Gavrila distribution for the K-shell is used to sample the
 * polar angle of a photoelectron. This performs the same sampling routine as
 * in Geant4's G4SauterGavrilaAngularDistribution class, as documented in
 * section 6.3.2 of the Geant4 Physics Reference (release 10.6) and section
 * 2.1.1.1 of the Penelope 2014 manual.
 */
template<class Engine>
CELER_FUNCTION Real3 PhotoelectricInteractor::sample_direction(Engine& rng) const
{
    constexpr MevEnergy min_energy{1.e-6};
    constexpr MevEnergy max_energy{100.};
    real_type           energy_per_mecsq;

    // If the incident gamma energy is above 100 MeV, use the incident gamma
    // direction for the direction of the emitted photoelectron.
    if (inc_energy_ > max_energy)
    {
        return inc_direction_;
    }
    // If the incident energy is below 1 keV, set it to 1 keV.
    else if (inc_energy_ < min_energy)
    {
        energy_per_mecsq = min_energy.value() * shared_.inv_electron_mass;
    }
    else
    {
        energy_per_mecsq = inc_energy_.value() * shared_.inv_electron_mass;
    }

    // Calculate Lorentz factors of the photoelectron
    real_type gamma = energy_per_mecsq + 1.;
    real_type beta  = std::sqrt(energy_per_mecsq * (gamma + 1.)) / gamma;
    real_type a     = (1. - beta) / beta;

    // Second term inside the brackets in Eq. 2.8 in the Penelope manual
    real_type b = 0.5 * beta * gamma * energy_per_mecsq * (gamma - 2.);

    // Maximum of the rejection function g(1 - cos \theta) given in Eq. 2.8,
    // which is attained when 1 - cos \theta) = 0
    real_type g_max = 2. * (1. / a + b);

    // Rejection loop: sample 1 - cos \theta
    real_type g;
    real_type one_minus_costheta;
    do
    {
        // Sample 1 - cos \theta from the distribution given in Eq. 2.9 using
        // the inverse function (Eq. 2.11)
        real_type u        = generate_canonical(rng);
        one_minus_costheta = 2. * a * (2. * u + (a + 2.) * std::sqrt(u))
                             / ((a + 2.) * (a + 2.) - 4. * u);

        // Calculate the rejection function (Eq 2.8) at the sampled value
        g = (2. - one_minus_costheta) * (1. / (a + one_minus_costheta) + b);
    } while (g < g_max * generate_canonical(rng));

    // Sample the azimuthal angle and calculate the direction of the
    // photoelectron
    UniformRealDistribution<real_type> sample_phi(0, 2 * constants::pi);
    return rotate(from_spherical(1. - one_minus_costheta, sample_phi(rng)),
                  inc_direction_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
