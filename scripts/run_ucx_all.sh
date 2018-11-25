#!/bin/bash

for ii in {1,2,4,8,16,28}
do
    for ss in {rc,rc_x,dc,dc_x}
    do
        for tt in {p,s,c,k}
        do
            if [[ $tt == k ]]
            then
                ll=100
                wl=10
            else
                ll=500
                wl=50
            fi
            ./run_ucx.sh -m=1:65536:2 -l=$ll -wl=$wl -b -type=$tt -ucx_tls=$ss -ppn=$ii
        done
    done
done
