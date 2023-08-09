#/bin/bash 
# Run this from the transfer folder. 

computer=${1:-brunsvik@Bellows.eri.ucsb.edu} # First argument is the computer you will send to. 

rsync -ahv \
--exclude-from='exclude_hk_anis.txt' \
../ \
$computer:~/Documents/repositories/hk_anis