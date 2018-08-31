#!/bin/bash

source activate cobra
TEMPLATE="./main_1.py"

BACKGROUND[0]="plastic"
BACKGROUND[1]="coating"
BACKGROUND[2]="czt"

if [ ! -f "$TEMPLATE" ]; then
 echo "Error: "$TEMPLATE" not found"
fi

for pyfiles in "${BACKGROUND[@]}"; do
  SKRIPT_NAME="main_1_"
  SKRIPT_NAME=$SKRIPT_NAME$pyfiles".py"
  FILE_NAME="main_1_"
  FILE_NAME=$FILE_NAME$pyfiles
    echo "creating script "$SKRIPT_NAME
    echo "creating file "$FILE_NAME
    if ! sed -e 's/eventfile_template/'eventfile_$pyfiles'/g' -e 's/scale_template/'scale_$pyfiles'/g' -e 's/x_range_template/'x_range_$pyfiles'/g' "$TEMPLATE" > "$SKRIPT_NAME"; then
      echo "Error: could not create Template"
    fi
  python $SKRIPT_NAME
done
python "./main_3.py"
source deactivate cobra
