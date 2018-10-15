#!/bin/bash

root -b -l <<-EOF
.L ./read_rootfile_energy_entries.cpp++
.L ./read_rootfile_energy_entries.cpp 
ausfuehren()
.q
EOF
