# DO NOT RUN ME - RUN ~/Git/UPSCb/src/bash/updateTaxonomySqlite.sh instead

#!/bin/env bash
set -ex

# go to the dir
cd /mnt/picea/storage/reference/Taxonomy

# create the dir
dat=`date "+%Y%m%d"`
mkdir $dat && cd $dat

# get the data
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz

# extract what we need
tar -zxf taxdump.tar.gz nodes.dmp names.dmp

# and reformart
sed -i 's:\t::g' nodes.dmp
sed -i 's:"::g' nodes.dmp
sed -i 's:|$::g' nodes.dmp
sed -i 's:\t::g' names.dmp
sed -i 's:"::g' names.dmp
sed -i 's:|$::g' names.dmp

# extract the rest
gunzip gi_taxid_nucl.dmp.gz
gunzip gi_taxid_prot.dmp.gz

# load the data in SQL
sqlite3 taxonomy.sqlite < $UPSCb/src/sql/taxonomy-create-and-populate-table.sql

# run R script to recreate the division
Rscript $UPSCb/src/R/updateTaxonomyDivisionTable.R 

# load the GI
sqlite3 taxonomy.sqlite < $UPSCb/src/sql/taxonomy-add-gi-relationship.sql

# compress
gzip gi_taxid_nucl.dmp gi_taxid_prot.dmp nodes.dmp names.dmp
