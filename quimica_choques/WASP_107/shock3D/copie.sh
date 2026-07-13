#!/bin/bash

for d in phi*; do
    cp -a shock_fiducial/. "$d"/
done
