#!/bin/bash

make clean
sphinx-apidoc -fMeET -o api-gen ../Inelastica
make html
