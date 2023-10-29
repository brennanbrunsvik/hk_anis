# Make HK stacks that incorporate radial anisotropy

To see this in action, download the package and run "example_workflow_for_github.m". 

Take seismic receiver functions and make HK stacks, while accounting for radial anisotropy. This package accompanies the paper Brunsvik and Eilon 2023, SRL, "Radial anisotropy in receiver function H−κ stacks". Contact me (brennanbrunsvik@ucsb.edu) if you need any help or would like my collaboration. 

- `example_workflow_for_github.m`: This is your starting point to understand the repository. Run it. Make both anisotropic and isotropic receiver functions for station TA KMSC. Example_workflow_complex.m has some more complicated options to handle more stations than just TA KMSC. 
- `hk_error_ray_param.m`: Calculate amount of error for a receiver function as a function of ray parameter. 
- hk_ears_bootstrap.m: Do bootstrapping to test uncertainty of HK stacks with station TA KMCS. 
- `summarize_many_sta_results.m`: If you were to run this on many stations, this script can be your starting place to plot many of those results in map view. 
- `hk_calculation/` has various files I used for calculating HK stacks. 
- `elastic/` has files used to calculate things like velocities of the model. 

To apply `example_workflow_for_github.m` to your own data, look for the comments "Load your data here!" and "Done loading your data!" as they indicate the variables you need to provide. You will need the following variables where nrf is the number of receiver functions, nt is the number of time samples, and nz is the number of depth points over which you are providing elastic parameters such as velocity. 
- `z`      [nz × 1 double]    (km) : Depth. 
  - I used values between 0 and 300. Values beneath zmoh are effectively ignored. 
- `vp`     [nz × 1 double]   (km/s) : P-wave velocity at each depth.
- `vs`     [nz × 1 double]   (km/s) : S-wave velocity at each depth.
- `xi`     [nz × 1 double]  : S-wave anisotropy at each depth.
  - Values will probably be between 0.8 and 1.2.
- `phi`    [nz × 1 double]  : P-wave anisotropy at each depth.
  - We used phi=1/xi in our example.
- `eta`    [nz × 1 double]  : Anisotropy parameter.
  - We used 1 for all depths. 
- `zmoh`   [scalar, double]   (km) : A first guess at the depth of the Moho.
  - We used 30 km. This value is used for applying the Earth flattening transformation to velocities down to this depth. It isn't needed if you already know your average crustal velocities and anisotropy. 
- `rf`     [nt × nrf double] : Your receiver functions.
  - You can choose your own style of normalization for the receiver functions, but the plotting code was built with the receiver functions having values between about -0.5 and 1.5. 
- `tt`     [nt × 1 double]   (seconds) : Time corresponding to each receiver function sample.
  - Each receiver function should share a common time basis.
- `rayp`   [1 × nt double]   (seconds/degree) : Ray parameter associated with each receiver function.
  - Units are in seconds per degree, with values typically between 4.5 and 8.8.

If you already have your crustal averaged elastic parameters, you can load your scalar values into the following variables and skip running the function `hk_average_crust`: `rhoav`, `vsav`, `vpav`, `xiav`, `phiav`, `etaav`. 

Synthetic receiver functions are not part of this package. It would be too much work to get our synthetic forward modeling working on other people's computers. You need to use the example receiver function included here, or provide your own. 
