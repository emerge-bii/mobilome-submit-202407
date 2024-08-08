#!/bin/bash

for f in fig/*.R; do
  echo "##############################"
  echo "[INFO] Processing $f"
  echo "##############################"
  Rscript $f || \
	  (echo && echo && \
	  echo "[ERROR] Failed: $f")

  echo
  echo
done
