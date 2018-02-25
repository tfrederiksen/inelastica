#!/bin/bash

make clean
sphinx-apidoc -fMeET -o api-gen ../Inelastica ../Inelastica/**/setup.py
make html
