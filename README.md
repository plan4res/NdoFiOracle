# The NDOSolver / FiOracle Project {#mainpage}


## Introduction

This is a splash page of the documentation of the `NDOSolver / FiOracle`
Project, a suite of ``C++`` interface classes and solvers for
NonDifferentiable Optimization problems (NDO) that have been developed at
the Department of Computer Science of the University of Pisa since about
1992 (although the current form of the interface, using the pair of
abstract classes NDOSolver and FiOracle, was only developed in 2001).

The structure of the software centers on a very clean separation between the
problem to be solved, i.e., a (possibly, sum-structured) nondifferentiable
function (possibly, with "easy" components) and simple (box) constraints,
and the algorithm that solves it. However, for specific applications like
[Lagrangian relaxation of structured
programs](http://pages.di.unipi.it/frangio/abstracts.html#AOR05),
information has to flow between the two not necessarily in the "obvious"
direction (variable values from solver to oracle, function results from
oracle to solver), but also in non-obvious ways.

The interface set by the `NDOSolver / FiOracle` pair allows for all the
necessary information exchange, starting with the  computation of
["convexified" primal
solutions](http://pages.di.unipi.it/frangio/abstracts.html#AOR05).
These [solve the original problem if it is
continuous](http://pages.di.unipi.it/frangio/abstracts.html#JOC99),
possibly providing [stronger
bounds](http://pages.di.unipi.it/frangio/abstracts.html#DAM01) and/or
[faster solution
times](http://pages.di.unipi.it/frangio/abstracts.html#MP05a) w.r.t.
standard methodologies, or can be used to [construct feasible
solutions](http://pages.di.unipi.it/frangio/abstracts.html#TPS03) for
hard, structured integer programs. This brings Lagrangian techniques to
parity with more established ones, like Linear Programming, in terms of
generated information, allowing to use standard [branching
rules](http://pages.di.unipi.it/frangio/abstracts.html#AOR05) and
separation routines for valid inequalities, which ultimately result
in [dynamic Lagrangian variables
generation](http://pages.di.unipi.it/frangio/abstracts.html#MP05a).

Yet, the interface cleanly separates the NDO solver from the function
(oracle), allowing to seamlessly test on the same problem algorithms
belonging to different classes, like different classes of [(generalized)
Bundle methods](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT02),
possibly [using nonstandard
models](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT11), and
[SubGradient-type
methods](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT08).
This is important because [different applications benefit from different
solvera](http://pages.di.unipi.it/frangio/abstracts.html#MPC16), as
well as on fine details about how the information gathered while computing
the function is used to construct its *model* that drive the solution
process. In particular, [more sophisticated
models](http://pages.di.unipi.it/frangio/abstracts.html#MP11b) can
significantly improve the convergence rate of algorithms, but may
incur in a significantly increased cost for the *master problem*. When the
cost of computing the solution is not very high, either naturally or
because one exploits the ability of several NDO algorithms to [only do that
approximately](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT16),
the trade-off between master problem time, convergence speed (number of
iterations) and function evaluation time becomes complex to navigate, and
it may also depend on the availability of [specialized methods to solve the
master problem](http://pages.di.unipi.it/frangio/abstracts.html#COR96).

The flexibility provided by the NDOSolver / FiOracle interface allows to
navigate these complex trade-offs by developing the problem-specific code
(the FiOracle) only once, and then deploying whatever general-purpose NDO
approach (the NDOSolver) is better suited for this task. The project
also provides state-of-the-art implementations of both main classes of
NDO algorithms, i.e., Bundle-type methods and SubGradient-type methods.


## NDO Solvers

Besides the abstract base classes `NDOSolver` and `FiOracle`, the project
provides the following algorithms:

- The `SubGrad` `NDOSolver` implements a large and parametric class of
  SubGradient-type algorithms that allow [simultaneous deflection and
  projection](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT08).
  The solver currently implements quite a large number of different stepsize
  and deflection rules, ranging from ad-hoc ones developed from Lagrangian
  relaxation of integer programs, to ones reproducing the [Volume
  algorithm](https://link.springer.com/article/10.1007/s101070050002), to
  full [primal-dual
  ones](https://link.springer.com/article/10.1007/s10107-007-0149-x).
  Although subgradient-type algorithms typically show slower convergence
  than cutting-plane based ones, their very low master problem cost per
  iteration still makes them [suitable for certain
  applications](http://pages.di.unipi.it/frangio/abstracts.html#MPC16),
  possibly in conjunction with [smoothing
  techniques](http://pages.di.unipi.it/frangio/abstracts.html#ORL17a).

- The `CutPlane` `NDOSolver` is a "didactic" implementation of [Kelley’s
  Cutting Plane approach](https://epubs.siam.org/doi/abs/10.1137/0108053),
  the ["father" of most modern NDO
  approaches](http://pages.di.unipi.it/frangio/abstracts.html#NDOB18).
  It relies on an [OsiSolverInterface](https://projects.coin-or.org/Osi)
  object to solve the master problem (a Linear Program), which allows the
  code to work with a large variety of both open-source and commercial
  solvers. It is provided more as an introduction to cutting-plane methods
  than as a strongly competitive algorithm, although in a few cases the
  "pure" Cutting Plane approach can actually turn out to be a reasonable
  choice due to the simpler master problem.

- The `Bundle` `NDOSolver` is a full-fledged implementation of a
  [generalized Bundle
  method](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT02),
  allowing for both proximal-type and trust-region-type stabilization terms.
  Indeed, it relies on the abstract class `MPSolver` to solve the master
  problem, which allows to abstract away both from the specific stabilization
  term employed and from the details of the solver used. The `Bundle` solver
  efficiently implements dynamic variables generation, re-optimization after
  changes in the Lagrangian subproblem, several schemes for allowing
  [inexact computation of the objective
  function](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT16), and
  supports the ["easy component"
  variant](http://pages.di.unipi.it/frangio/abstracts.html#MP11c). It is
  currently distributed with two specific `MPSolver` implementations:

  = `QPPenaltyMP`, a [custom-made active-set
    solver](http://pages.di.unipi.it/frangio/abstracts.html#COR96) for the
    master problem of the Proximal Bundle method with one component;

  = `OSIMPSolver`, which relies on an
    [OsiSolverInterface](https://projects.coin-or.org/Osi) object to solve
    the master problem with either proximal or trust region stabilization,
    supporting disaggregated models for FiOracle with multiple components
	as well as ["easy
	components"](http://pages.di.unipi.it/frangio/abstracts.html#MP11c).

  From the user's perspective, selecting the `MPSolver` amounts at just
  changing a very few lines of code in the initialization of the `Bundle`
  solver, and possibly tuning a few corresponding algorithmic parameters.
  Yet, the two MPSolver provide [distinct advantages in different
  applications](http://pages.di.unipi.it/frangio/abstracts.html#MP11c),
  which significantly improve the capability of the Bundle solver to
  efficiently adapt to the requirements of the specific problem to be
  solved.


## Test FiOracle

A few very simple implementation of "test" `FiOracle` are provided as a
guide for developers of actual `FiOracle`. Other examples of `FiOracle`,
in particular related to the Fixed-Charge Multicommodity Capacitated
Network Design problem [are publicly
available](https://github.com/frangio68/SubGrad-4-FC-MMCF), and several
others still have been developed over the years and can be requested to
the developers of the code.

- `TestFi` implements an exceedingly simple convex quadratic function
  with one component;

- `LukFi` implements a single `FiOracle` interface for a set of
  [25 "standard" problems for nonsmmoth unconstrained
  optimization](http://www.cs.cas.cz/~luksan/test.html), which can be
  useful for benchmarking NDO solvers.


## Documentation

Full Doxygen documentation can be produced by just running doxygen into
the [`doxygen/`](doxygen) folder. It is [also available online
here](https://frangio68.gitlab.io/ndosolver_fioracle_project).


## Software Dependencies

The `NDOSolver / FiOracle` interface is written in standard (pre-11) ``C++``,
with no specific requirements for external libraries. The `Bundle` solver with
the `QPPenaltyMP` master problem solver and the `SubGrad` solver are similarly
completely unencumbered. The `CutPlane` solver and the `OSIPSolver` master problem
solver depend on the [OsiSolverInterface](https://projects.coin-or.org/Osi)
project from COIN-OR; they can use the free
[Clp](https://projects.coin-or.org/Clp) solver from COIN-OR or almost any
commercial solver like Cplex, GuRoBi, Mosek, ... (to be separately obtained
and licensed).


## Installation Instructions

To download a copy of the software, just run

```sh
git clone https://gitlab.com/frangio68/ndosolver_fioracle_project
```

You can either use [CMake](https://cmake.org) or plain makefiles to build the
library, your choice. CMake compiles off-source and it is therefore perhaps
better suited to one-off, compile-and-forget installations, whereby the
provided makefiles compile on-source and we find that they are better suited
while developing and testing the code (if that's your cup of tea; it is ours).

In both cases, all external dependencies should be automatically dealt with if
they are installed in their default paths, as specified in the `*_ROOT` values
of [`extlib/makefile-default-paths`](extlib/makefile-default-paths). If not,
the suggested way to change them is to copy the file into
[`extlib/makefile-paths`](extlib/makefile-paths) and edit it. The file (if
present) is automatically read and the values found there replace the
corresponding non-default definitions. The rationale for not changing
makefile-default-paths is that makefile-paths file is .gitignore-d. Hence, it
should not be necessary to re-change the makefiles (or stash/restore the
changes) each time the project is pulled, or manually ignore the changes when
it is pushed, which is very convenient for anyone who actually develops
`NDOSolver / FiOracle` components (anyone there?). However, note that the
reading of both
[`extlib/makefile-default-paths`](extlib/makefile-default-paths) and
[`extlib/makefile-paths`](extlib/makefile-paths) can be disabled; see the
`NDOFiOracle_READ_PATHS` option in the CMake section and the `NDOFi_NO_PATHS`
macro in the makefile section below.

### Using CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

If the required libraries are not at their default location you can specify
actual install location via the
[`extlib/makefile-paths`](extlib/makefile-paths) file, as specified above.

You can also choose the following configuration options:

| Variable              | Description         | Default value |
|-----------------------|---------------------|---------------|
| `NDOFiOracle_USE_CLP` | Use Clp             | ON            |
| `NDOFiOracle_USE_OSI` | Use Osi             | ON            |
| `WHICH_OSI_MP`        | Use CPLEX or GUROBI | GUROBI        |

You can set them with:

```sh
cmake <source-path> -D<var>=<value>
```

Optionally, install the library in the system with:

```sh
sudo make install
```

After the library is built, you can use it in your CMake project with:

```cmake
find_package(NDOFiOracle)
target_link_libraries(<my_target> NDOFiOracle)
```

- If you use the `NDOSolver / FiOracle` project inside some other project that
  already properly defines the `*_ROOT` values, you can avoid them being read by
  setting the option `NDOFiOracle_READ_PATHS` to `OFF`, e.g., by adding

  ```cmake
  set(NDOFiOracle_READ_PATHS OFF CACHE BOOL
                          "Whether NDOSolver / FiOracle will read locations for
                          dependencies or not." FORCE)
  ```

  in your CMake project.

### Installation Instructions for makefiles

The main makefile is [`lib/makefile`](lib/makefile). There, one finds the main
sub-modules of the project. The `NDOSlver` sub-module is mandatory, as it defines
the `NDOSolver / FiOracle` interface. All other sub-modules define specific
algorithmic approaches and are optional. In particular:

- `Bundle` is the generalized-bundle-type algorithm. It relies on separate
  sub-modules (implementing the `MPSolver` class) to solve the "master problem":

  = `QPPnltMP`, using a taylor-made QP solver

  = `OSIMPSolver`, using a generic `OsiSolverInterface`

  If the `Bundle` sub-module is selected (not commented out), at least one of
  these must be selected. While `QPPnltMP` is fully self-contained, `OSIMPSolver`
  requires a working `OsiSolverInterface` distribution, as detailed below.
  Also, the `MPTester` sub-module is provided that can be useful for developers
  of alternative `MPSolver` since it simplifies comparing results of
  different `MPSolver`.

- The `CutPlane` sub-module is a "didactic" implementation of the Cutting Plane
  approach, probably not much useful in most cases. It also requires a working
  OsiSolverInterface distribution.

- The `SubGrad` sub-module which implements a number of variants of subgradient
  approaches, which is fully self-contained.

Once you have selected which modules you want to be compiled in the library,
by just un-commenting the corresponding include lines, if these are not fully
self-contained you will have to be sure that the paths for the dependencies
are properly set. This could be free if they are installed at their default
locations, otherwise you can use the
[`extlib/makefile-paths`](extlib/makefile-paths) file, as specified above. Note
that you also have to be ensure that the corresponding include line in
[`lib/makefile-c`](lib/makefile-c) is uncommented.

However, a specific twist also need to be taken care of for `OSIMPSolver`. This
is because the `OsiSolverInterface` does not support solving quadratic programs,
while many underlying solvers actually allow it, and this is useful for the
bundle approach (since quadratic stabilizations are often preferable to linear
ones). Thus, the compile-time switch

```C
#define WHICH_OSI_MP
```

in `OSIMPSolver.C` activates solver-specific code that allows to use quadratic
stabilization with selected solvers, in particular Cplex and GuRoBi. However,
you then need

- to have the corresponding `OsiCpxSolverInterface` or `OsiGrbSolverInterface`
  modules built in the `OsiSolverInterface` (cf. the documentation of that
  project for details)

- to ensure that the corresponding paths are added while building the
  library, which is done similarly to `OsiSolverInterface` by properly
  editing [`lib/makefile-c`](lib/makefile-c) and ensuring that the
  corresponding include line in `lib/makefile-c` is not commented out

Once this is done, building the library just (possibly) requires to edit
[`lib/makefile-lib`](lib/makefile-lib) for little details like the compiler
name and the optimization/debug switches, and then run

```sh
cd lib
make -f makefile-lib 
```

This creates `lib/libNDO.a` that can be linked upon. However, note that the
include files still remain "scattered" across the source directories, as the
makefile does not provide an "install" option to gather them in a unique place.
However, the file [`lib/makefile-inc`](lib/makefile-inc) is provided for this
purpose. The simplest way to learn how to use it is to check the two examples
of testing executables, located in [`LukFi/test`](LukFi/test) and
[`TestFi/test`](TestFi/test). The executables can be built by simply running
"make" there, which also builds the library. Note that, in this case, compiler
name and options are rather chosen in the "main" makefile of the test
executable. The FiOracles, main.C and makefiles in the two examples also provide
good starting points for developing your own `FiOracle`, and supporting
makefiles/test executables to use it.


## Authors

### Current Lead Author

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Enrico Gorgone**  
  Dipartimento di Matematica e Informatica  
  Università di Cagliari

- **Vitor Barbosa**  
  Algoritmi Research Unit  
  Universidade do Minho

- **Filipe Alvelos**  
  Departamento de producao e Sistemas  
  Universidade do Minho

- **Luigi Poderico**  

- **Andrea Nerli**  

- **Annabella Astorino**  
  Istituto di Calcolo e Reti ad Alte Prestazioni 
  Consiglio Nazionale delle Ricerche

- **Manlio Gaudioso**  
  Dipartimento di Elettronica Informatica e Sistemistica  
  Università della Calabria


## License

This code is provided free of charge under the "GNU Lesser General Public
License" version 3.
