#!/bin/bash

DATA_FILE_PATH=$1
SPECIES=$2
echo "Input file: $DATA_FILE_PATH"
echo "Species: $SPECIES"

gunzip -c $DATA_FILE_PATH | grep "OS=$SPECIES" \
| cut -d '|' -f 2 | wc -l
