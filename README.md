# GradientSensing

Code base for:
*Slower swimming promotes chemotactic encounters between bacteria and small phytoplankton*,
R. Foffi, D.B. Brumley, F.J. Peaudecerf, R. Stocker, J. Słomka.

The code base is authored by Riccardo Foffi.

## How to reproduce
To (locally) reproduce this project, do the following:

0. Install Julia () and download this code base. Notice that no data is included in the git-history.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Data production pipeline
- Run `scripts/parameterspace_hein.jl` to produce the range of $R$ and $C_S$ values that
  will be used in the numerical evaluation of the theoretical $I_C$ landscape.
  It will output `RC.jld2` to `datadir("Hein")`
- Run `scripts/phycosphere_hein_mod.jl` to numerically evaluate the sensory horizon $S$ for the selected set of parameters.
  For each parameter combination it will output a file with prefix `phycosphere` to `datadir("Hein")`.
- Run `scripts/ic_hein.jl` to evaluate the $I_C$ using the numerical computations of $S$.
  For each parameter combination it will output a file with prefix `IC` to datadir("Hein").
- Run `scripts/parameterspace_poisson.jl` to produce the range of $R$ and $C_S$ values that
  will be used in the evaluation of the $I_C$ landscape through the Kolmogorov-Smirnov sensor.
  It will output `RC.jld2` to `datadir("Poisson")`
- Run `scripts/poisson_sampling.jl` to simulate spatial sampling of concentration profiles
  by measuring waiting times between adsorption events along a transect towards a spherical
  source.
  For each parameter combination it will output a file with prefix `waitingtimes` to `datadir("Poisson")`.
- Run `scripts/sensing_kolmogorovsmirnov.jl` to perform Kolmogorov-Smirnov tests of the distributions of
  waiting times sampled by `poisson_sampling.jl`.
  For each parameter combination it will output a file with prefix `sensing` to `datadir("Poisson", "KolmogorovSmirnov")`.
- Run `scripts/phycosphere_kolmogorovsmirnov.jl` to evaluate the sensory horizon $S$ associated with the
  Kolmogorov-Smirnov tests.
  For each parameter combination it will output a file with prefix `phycosphere` to `datadir("Poisson", "KolmogorovSmirnov")`.
- Run `scripts/ic_kolmogorovsmirnov.jl` to evaluate the $I_C$ associated with the
  Kolmogorov-Smirnov estimates of the sensory horizon $S$.
  For each parameter combination it will output a file with prefix `IC` to `datadir("Poisson", "KolmogorovSmirnov")`.

## Dashboards
The dashboards for the interactive evaluation of the $I_C$ landscape and of
bacteria-phytoplankton search times are found in the `dashboards` directory.

It is not required to run the data production pipeline to use the dashboards.
