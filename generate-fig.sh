#!/bin/bash

for f in fig/*.R; do
  Rscript $f
done
