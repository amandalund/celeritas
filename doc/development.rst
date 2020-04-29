=====================
Celeritas development
=====================


Core guidelines
===============

In which the author writes his first manifesto since probably high school.

Maximize encapsulation
----------------------

Encapsulation is about making a piece of code into a black box. The fewer lines
connecting these black boxes, the more maintainable the code. Black boxes can
often be improved internally by making tiny black boxes inside the larger black
box.

Motivation:

- Developers don't have to understand implementation details when looking at a
  class interface.
- Compilers can optimize better when dealing with more localized components.
- Good encapsulation allows components to be interchanged easily because they
  have well-defined interfaces.
- Pausing to think about how to minimize input and output from an algorithm can
  improve make it easier to write.

Applications:

- Refactor large functions (> 50 statements ish?) into small functors that take
  "invariant" values (the larger context) for constructors and use
  ``operator()`` to transform some input into the desired output
- Use only ``const`` data when sharing. Non-const shared data is almost like
  using global variables.
- Use ``OpaqueId`` instead of integers and magic sentinel values for
  integer identifiers that aren't supposed to be arithmetical .

Examples:

- Random number sampling: write a unit sphere sampling functor instead of
  replicating a polar-to-cartesian transform in a thousand places
- Cell IDs: Opaque IDs add type safety so that you can't accidentally convert a
  cell identifier into a double or switch a cell and material ID. Also makes
  code more readable of course.


Maximize code reuse
-------------------

No explanation needed.


Minimize compile time
---------------------

Code performance is important, but so is developer time. When possible,
minimize the amount of code touched by NVCC. (NVCC's error output is also
rudimentary compared to modern clang/gcc, so that's another reason to prefer
them compiling your code.)

Prefer single-state classes
---------------------------

As much as possible, make classes "complete" and valid after calling the
constructor. Don't have a series of functions that have to be called in a
specific order to put the class in a workable state.

When a class has a single function (especially if you name that function
``operator()``), its usage is obvious. The reader also doesn't have to know
whether a class uses ``doIt`` or ``do_it`` or ``build``.

When you have a class that needs a lot of data to start in a valid state, use a
``struct`` of intuitive objects to pass the data to the class's constructor.
The constructor can do any necessary validation on the input data.


Style guidelines
================

Having a consistent code style makes it more readable and maintainable. You
don't have to guess whether a function or class.

File names
----------

All ``__device__`` and ``__global__`` code must be compiled with NVCC to generate
device objects. However, code that merely uses CUDA API calls such as
``cudaMalloc`` does *not* have to be compiled with NVCC. Instead, it only has to
be linked against the CUDA runtime library and include ``cuda_runtime_api.h``.
The exception to this is VecGeom's code, which compiles differently when run
through NVCC. (Macro magic puts much of the code in a different namespace.)

Since NVCC is slower and other compilers' warning/error output is more
readable, it's preferable to use NVCC for as little compilation as possible.
Furthermore, not requiring NVCC lets us play nicer with downstream libraries
and front-end apps. Host code will not be restricted to the minimum version
supported by NVCC (C++14).

Of course, the standard compilers cannot include any CUDA code containing
kernel launches, since those require special parsing by the compiler. So kernel
launches and ``__global__`` code must be in a ``.cu`` file. However, the
CUDA runtime does define the special ``__host__`` and ``__device__`` macros (among
others). Therefore it is OK for a CUDA file to be included by host code as long
as it ``#include`` s the CUDA API. (Note that if such a file is to be included by
downstream code, it will also have to propagate the CUDA include directories.)

Choosing to compile code with the host compiler rather than NVCC also provides
a check against surprise kernel launches. For example, the declaration::

   thrust::device_vector<double> dv(10);

actually launches a kernel to fill the vector's initial state. The code will
not compile in a ``.cc`` file run through the host compiler, but it will
automatically (and silently) generate kernel code when run through NVCC.

Thus we have the following rules:

- ``.h`` is for C++ code compatible with host compilers. It may declare Thrust
  objects, since thrust type declarations are compatible with the host
  compiler. It can also use host/device keywords if it includes the cuda
  runtime api or hides the keywords with macros.
- ``.cu`` is for ``__global__`` kernels and functions that launch them
- ``.cuh`` is for header files that require compilation by NVCC: contain
  ``__device __``-only code or include CUDA directives without ``#include
  <cuda_runtime_api.h>``.

Code object names
-----------------

Functions should be verbs; classes should be names. Functors (classes whose
instances act like a function) should be an *agent noun*: the noun form of an
action verb. Instances of a functor should be a verb. For example::

   ModelEvaluator evaluate_something(parameters...);
   auto result = evaluate_something(arguments...);

Variable names
--------------

Generally speaking, variables should have short lifetimes and should be
self-documenting. Avoid shorthand and "transliterated" mathematical
expressions: prefer ``constants::avogadro`` to ``N_A`` or express the constant
functionally with ``atoms_per_mole``.


