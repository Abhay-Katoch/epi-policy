#!/bin/bash

directory="/Users/abhay/Documents/XLab/epi-policy"

echo "Running results.py"
python "${directory}/epi-policy/results.py"

echo "Running figures.py"
python "${directory}/epi-policy/figures.py"