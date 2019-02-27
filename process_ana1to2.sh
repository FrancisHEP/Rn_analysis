#!/bin/zsh
# ls /nas/64-core/$1/ana_run$1_pattern_file* > filelist/ana_run$1.lst
# ls /nas/64-core/$1/run$1_pattern_file*hit.root > filelist/hit_run$1.lst
# /home/andy/pandax-chain/build/ana1to2 filelist/hit_run$1.lst filelist/ana_run$1.lst ana2/run$1.root
run=$1
ls /store/px/data/udm/by-run/${run}/new_correction/ana_run${run}_pattern_file* > filelist/ana_run${run}.lst
ls /store/px/data/udm/by-run/${run}/new_correction/run${run}_pattern_file*hit.root > filelist/hit_run${run}.lst
/home/mawenbo/pandax-chain/build/ana1to2 filelist/hit_run${run}.lst filelist/ana_run${run}.lst ana2/run${run}.root

