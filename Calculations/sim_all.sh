PRAEFIX[0]="CZT_double_beta_decay_"
PRAEFIX[1]="Coatin_chain_decay_"
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
  
done

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
