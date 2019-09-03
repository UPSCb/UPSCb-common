#!/bin/bash
set -e

# set the dir
DIR=/mnt/picea/storage/reference/Taxonomy/`date "+%Y%m%d"`
mkdir $DIR
cd $DIR

# retrieve the data
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz

# extract what we need
tar -zxf taxdump.tar.gz nodes.dmp names.dmp
find . -name "gi*.dmp.gz" | xargs -P 2 -I{} gunzip {}

# update some table dumps
sed -i 's:\t::g' nodes.dmp
sed -i 's:"::g' nodes.dmp
sed -i 's:|$::g' nodes.dmp
sed -i 's:\t::g' names.dmp
sed -i 's:"::g' names.dmp
sed -i 's:|$::g' names.dmp

# create and populate the database
sqlite3 taxonomy.sqlite < $UPSCb/src/sql/taxonomy-update.sql

# update the database dynamically (fix some awkwardities in the taxonomy tables)
module load R
Rscript ~/Git/UPSCb/src/R/updateTaxonomyDivisionTable.R

# clean up
rm names.dmp nodes.dmp
find . -name "*.dmp" | xargs -P 2 -I{} gzip {}
