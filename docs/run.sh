#!/bin/bash

rm -rf api-gen
make clean
sphinx-apidoc -fMeET -o api-gen ../Inelastica ../Inelastica/**/setup.py
make html
