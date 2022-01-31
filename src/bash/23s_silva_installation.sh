#!/bin/bash

# Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build a 16S database from Silva data

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

HTTPS_SERVER="https://ftp.arb-silva.de/"
SILVA_VERSION="138_1"
REMOTE_DIR="$HTTPS_SERVER/release_$SILVA_VERSION/Exports"
FASTA_FILENAME="SILVA_${SILVA_VERSION/_/.}_LSURef_NR99_tax_silva.fasta"
TAXO_PREFIX="tax_slv_lsu_${SILVA_VERSION/_/.}"

if [ $# -ne 1 ]; then
  echo "Set the output path as the first argument"
  exit 1
fi

export PATH=$PATH:$(realpath ../../kogia/scripts)
#export PATH=$PATH:$(realpath ../../../kraken2/scripts)

KRAKEN2_DB_NAME=$1
KRAKEN2_THREAD_CT=1


mkdir -p "$KRAKEN2_DB_NAME"
#pushd "$KRAKEN2_DB_NAME"
cd $KRAKEN2_DB_NAME
mkdir -p data taxonomy library
#pushd data
cd data
if [ ! -f ${FASTA_FILENAME} ]; then
  wget -nc "$REMOTE_DIR/${FASTA_FILENAME}.gz"
  gunzip -f "${FASTA_FILENAME}.gz"
fi
if [ ! -f ${TAXO_PREFIX}.acc_taxid ]; then
  wget -nc "$REMOTE_DIR/taxonomy/${TAXO_PREFIX}.acc_taxid.gz"
  gunzip -f "${TAXO_PREFIX}.acc_taxid.gz"
fi
if [ ! -f ${TAXO_PREFIX}.txt ]; then
  wget -nc "$REMOTE_DIR/taxonomy/${TAXO_PREFIX}.txt.gz"
  gunzip -f "${TAXO_PREFIX}.txt.gz"
fi

#docker run --rm --entrypoint /usr/local/kraken2/build_silva_taxonomy.pl -v /mnt:/mnt bschiffthaler/kraken2 $KRAKEN2_DB_NAME/data/"${TAXO_PREFIX}.txt"
if [ ! -f $KRAKEN2_DB_NAME/taxonomy/names.dmp ] || [ ! -f $KRAKEN2_DB_NAME/taxonomy/nodes.dmp ]; then
  $(realpath ../../../kraken2/scripts/build_silva_taxonomy.pl) "${TAXO_PREFIX}.txt"

#popd
  mv $KRAKEN2_DB_NAME/data/names.dmp $KRAKEN2_DB_NAME/data/nodes.dmp $KRAKEN2_DB_NAME/taxonomy/
  mv $KRAKEN2_DB_NAME/data/${TAXO_PREFIX}.acc_taxid $KRAKEN2_DB_NAME/seqid2taxid.map
fi

if [ ! -f $KRAKEN2_DB_NAME/library/silva.fna ]; then
  sed -e '/^>/!y/U/T/' "$KRAKEN2_DB_NAME/data/$FASTA_FILENAME" > $KRAKEN2_DB_NAME/library/silva.fna
fi
#popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT
