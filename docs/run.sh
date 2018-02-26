#!/bin/bash

rm -rf api-gen latest
make clean
sphinx-apidoc -fMeET -o api-gen ../Inelastica ../Inelastica/**/setup.py
make html
rm -r build/html/_sources
mv build/html latest
