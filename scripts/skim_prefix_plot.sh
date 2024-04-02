#!/bin/bash

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <file.blow5>"
    exit 1
fi

slow5tools skim --rid ${1} | while read p 
do
    ./scripts/plot_prefix.sh ${1} ${p}
done