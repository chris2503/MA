#!/bin/bash

source activate cobra
TEMPLATE="./compare_versions.py"

BACKGROUND[0]="plastic"
BACKGROUND[1]="coating"
BACKGROUND[2]="czt"

if [ ! -f "$TEMPLATE" ]; then
 echo "Error: "$TEMPLATE" not found"
fi

for pyfiles in "${BACKGROUND[@]}"; do
  SKRIPT_NAME="compare_versions_"
  SKRIPT_NAME=$SKRIPT_NAME$pyfiles".py"
  FILE_NAME="compare_versions_"
  FILE_NAME=$FILE_NAME$pyfiles
    echo "creating script "$SKRIPT_NAME
    echo "creating file "$FILE_NAME
    if ! sed -e 's/eventfile_template_1/'eventfile_${pyfiles}_1'/g' -e 's/eventfile_template_2/'eventfile_${pyfiles}_2'/g' -e 's/scale_template/'scale_$pyfiles'/g' -e 's/x_range_template/'x_range_$pyfiles'/g' -e 's/background_template/'background_$pyfiles'/g' -e 's/background_2_template/'background_2_$pyfiles'/g' "$TEMPLATE" > "$SKRIPT_NAME"; then
      echo "Error: could not create Template"
    fi
  python $SKRIPT_NAME
done

source deactivate cobra
