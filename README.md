# Make HK stacks that incorporate radial anisotropy

This repository will be re-cleaned and organized upon acceptance of the corresponding manuscript. Right now, some script names may not be as expected. If you are trying to use this repository before acceptance of the manuscript, then after downloading the IRIS EARS data, you can run hk_ears_loop.m and show_many_stats.m to reproduce some figures from the paper. 

Take seismic receiver functions and make HK stacks, while accounting for radial anisotropy. This package accompanies the paper "Radial anisotropy in receiver function H−κ stacks" by Brennan Brunsvik and Zach Eilon, submitted to SRL on 2023/04/13. Please feel free to contact me (brennanbrunsvik@ucsb.edu) if you need any help or would like my collaboration. We would be very happy if somebody wanted to convert all of this to Python. 

- hk_ears.m: Make anisotropic HK stack on receiver functions from the IRIS-EARS dataset. Can loop over many stations, or just do one station. You can modify this to use receiver functions from any station. 
  - hk_ta_kmsc.m: Temporary, make anisotropic HK stack for a unique format station, data output was from our MCMC work. 
  - Note that we do not include the receiver function waveforms from IRIS-EARS. Instead, include the code to download it. For simple instructions, see models/download_rfs/README.txt. 

- hk_error_ray_param.m: Calculate amount of error for a receiver function as a function of ray parameter. 

Note that synthetic receiver functions are not part of this package. It would be too much work to get our synthetic forward modeling working on other people's computers, but if you have a way to calculate your own synthetic receiver functions, you can definitely use this package to make anisotropic HK stacks on them! 
