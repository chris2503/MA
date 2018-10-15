#!/bin/bash

source activate cobra

./analyse_p1.sh

python "./main_2.py"

python "./main_3.py"
source deactivate cobra
