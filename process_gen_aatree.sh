#!/bin/zsh
# ls /nas/64-core/$1/ana_run$1_pattern_file* > filelist/ana_run$1.lst
# ls /nas/64-core/$1/run$1_pattern_file*hit.root > filelist/hit_run$1.lst
# /home/andy/pandax-chain/build/ana1to2 filelist/hit_run$1.lst filelist/ana_run$1.lst ana2/run$1.root
run=$1
root -l -b -q 'Rn220_Po216.C('${run}')'

