# Make HK stacks that incorporate radial anisotropy

To see this in action, download the package and run "example_workflow_for_github.m". 

Take seismic receiver functions and make HK stacks, while accounting for radial anisotropy. This package accompanies the paper Brunsvik and Eilon 2023, SRL, "Radial anisotropy in receiver function H−κ stacks". Contact me (brennanbrunsvik@ucsb.edu) if you need any help or would like my collaboration. 

- example_workflow_for_github.m: This is your starting point to understand the repository. Run it. Make both anisotropic and isotropic receiver functions for station TA KMSC. Example_workflow_complex.m has some more complicated options to handle more stations than just TA KMSC. 
- hk_error_ray_param.m: Calculate amount of error for a receiver function as a function of ray parameter. 
- hk_ears_bootstrap.m: Do bootstrapping to test uncertainty of HK stacks with station TA KMCS. 
- summarize_many_sta_results.m: If you were to run this on many stations, this script can be your starting place to plot many of those results in map view. 
- hk_calculation/ has various files I used for calculating HK stacks. 
- elastic/ has files used to calculate things like velocities of the model. 

Synthetic receiver functions are not part of this package. It would be too much work to get our synthetic forward modeling working on other people's computers. You need to use the example receiver function included here, or provide your own. 
