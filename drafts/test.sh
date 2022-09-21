#!/bin/bash

input_files=(
    anna.col
    brock200_2.clq
    brock200_4.clq
    p_hat300-2.clq
)

if make -C build && rm -f max_clique && ln -s build/max_clique max_clique;
then
    for f in "${input_files[@]}";
    do
        echo -n "$f: "
        time ./max_clique < $f
        echo
    done
fi
