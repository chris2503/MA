#!/bin/bash

PRAEFIX[0]="CZT_double_beta_decay_"
PRAEFIX[1]="Coating_chain_decay_"
PRAEFIX[2]="Delrin_Layer_"
PRAEFIX[3]="Delrin_screws_"

CZT_ISOTOPES[0]="114Cd"
CZT_ISOTOPES[1]="116Cd"
CZT_ISOTOPES[2]="70Zn"
CZT_ISOTOPES[3]="128Te"
CZT_ISOTOPES[4]="130Te"

CHAIN_ISOTOPES[0]="232Th"
CHAIN_ISOTOPES[1]="238U"
CHAIN_ISOTOPES[2]="40K"

COATING[0]="_epoxy"
COATING[1]="_glyptal"

SUEFFIX=".mac"


for background in "${PRAEFIX[@]}"; do
  FILE_NAME=$background
  if [ "$background" == "${PRAEFIX[0]}" ]; then
    for isotopes in "${CZT_ISOTOPES[@]}"; do
      FILE_EXE="./"$FILE_NAME$isotopes$SUEFFIX
      echo "Executing venom "$FILE_EXE" ...
        "
      #parallel -j+0 --eta 
      venom $FILE_EXE > temp.txt
    done
  elif [ "$background" == "${PRAEFIX[1]}" ]; then
    for ch_isotopes in "${CHAIN_ISOTOPES[@]}"; do
      for coatings in "${COATING[@]}"; do
        FILE_EXE="./"$FILE_NAME$ch_isotopes$coatings$SUEFFIX
        echo "Executing venom "$FILE_EXE" ...
      "
        #parallel -j+0 --eta 
        venom $FILE_EXE > temp.txt
      done
    done
  else
    for ((COUNTER=0;COUNTER<2;COUNTER++)); do
      ch2_isotopes=${CHAIN_ISOTOPES[COUNTER]}
      FILE_EXE="./"$FILE_NAME$ch2_isotopes$SUEFFIX
      echo "Executing venom "$FILE_EXE" ...
      "
      #parallel -j+0 --eta 
      venom $FILE_EXE > temp.txt
    done
  fi
done

rm -f temp.txt
