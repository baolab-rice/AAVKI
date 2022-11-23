# AAVKI

touch AAV_Arti_Genome.fa

cat mm10.fa > AAV_Arti_Genome.fa
# AAV_CRISPR
cat AAV_CRISPR.fa >> AAV_Arti_Genome.fa 
# AAV_Donor
cat AAV_Donor.fa >> AAV_Arti_Genome.fa 
#ITR only to filter out virus recomb or NHEJ integration 
## Do we need to remove ITR from other template? NEED TO CHECK
##! Yes. Re-prepare the AAV templates w/o ITR
cat ITR.fa >> AAV_Arti_Genome.fa 
