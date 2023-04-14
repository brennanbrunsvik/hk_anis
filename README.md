# Make HK stacks that incorporate radial anisotropy
Note: As of 2023/04/13, this repository is made public. However, I need to clean it up, push a last commit with some edits, etc. If you are trying to access this and still see this message, please contact me. 

Take seismic receiver functions and make HK stacks, while accounting for radial anisotropy. This package accompanies the paper XXXXXXXXXXXXXXXXXXXXXX. Please feel free to contact me (brennanbrunsvik@ucsb.edu) if you need any help or would like my collaboration. We would be very happy if somebody wanted to convert all of this to Python, or a fast compiled language. 
 
% hk_figs_station_paper.m shows an examples of making anisotropic HK stacks for a particular station. It makes some figures used in the paper. 

% hk_error_paper.m makes some figures that show, analytically, how much error can be induced from ignoring anisotropy. 

hk_ears.m: Make anisotropic HK stack on receiver functions from the IRIS-EARS dataset. Can loop over many stations, or just do one station. 