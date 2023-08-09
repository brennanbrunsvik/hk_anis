#/bin/bash 
# Run this from the transfer folder. 

computer=${1:-brunsvik@Bellows.eri.ucsb.edu} # First argument is the computer you will send to. 

rsync -ahv \
/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/repositories/hk_anis/models/download_rfs/Ears/gauss_2.5/AE.X16A/** \
$computer:~/Documents/repositories/hk_anis/models/download_rfs/Ears/gauss_2.5/AE.X16A/