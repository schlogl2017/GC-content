#!/bin/bash

FNAME="data/uniprot_sprot.$(date +'%Y-%m-%d').fasta.gz"
curl --location --output $FNAME http://bit.ly/1l6SAKb

