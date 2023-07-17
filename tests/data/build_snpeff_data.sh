#!/usr/bin/env bash

# Use a conda environment with snpEff installed

snpEff build \
  -c snpEff.config \
  -dataDir snpEff_data \
  -v \
  -nodownload \
  testref1
