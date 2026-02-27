#!/usr/bin/env bash

# bash wrapper script for running SPfast

# WARN: be sure SPfast conda env is activated

# Precompute representative pseudoatoms
python ./utils/idealize.py -h

python ./utils/idealize.py \
    --structure_suffix ent.gz \
    --sdir ./example/structures \
    --odir ./output \
    --af2model \
    --trim \
    ./example/example_list # --dssdir ./example/DSSP \
