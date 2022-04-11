#!/bin/bash


DATABASE="photorhabdus_xenorhabdus_test.db"
REF_GENOMES_DIR="actinobacteria_internal streptomyces_internal"
QUERY_GENOMES_DIR="photorhabdus xenorhabdus"


### Fist part on local machine, as the query genomes are not on the shared machine
### to use the O drive locally (ubuntu subsystem on win) the drive needs to be mounted first: sudo mount -t drvfs O: /mnt/o

python -m cProfile -o ../results/profs/01_create_database_test.prof \
    01_create_database.py -db $DATABASE \
    2>&1 | tee ../results/logs/01_create_database_test.log 

python -m cProfile -o ../results/profs/02_import_genomes_ref_test.prof \
    02_import_genomes.py -db $DATABASE -r $REF_GENOMES_DIR --test \
    2>&1 | tee  ../results/logs/02_import_genomes_ref_test.log


### Second part on shared machine

python -m cProfile -o ../results/profs/02_import_genomes_query_test.prof \
    02_import_genomes.py -db $DATABASE -q $QUERY_GENOMES_DIR --test \
    2>&1 | tee  ../results/logs/02_import_genomes_query_test.log

python -m cProfile -o ../results/profs/03_add_core_genome_indicator_test.prof \
    03_add_core_genome_indicator.py -db $DATABASE \
    2>&1 | tee  ../results/logs/03_add_core_genome_indicator_test.log

python -m cProfile -o ../results/profs/04_find_approximate_pfamsequences_test.prof \
    04_find_approximate_pfamsequences.py -db $DATABASE -t 2 \
    2>&1 | tee  ../results/logs/04_find_approximate_pfamsequences_test.log

python -m cProfile -o ../results/profs/05_cluster_hits_test.prof \
    05_cluster_hits.py -db $DATABASE -t 2 \
    2>&1 | tee  ../results/logs/05_cluster_hits_test.log

python -m cProfile -o ../results/profs/06_add_cluster_information_test.prof \
    06_add_cluster_information.py -db $DATABASE \
    2>&1 | tee  ../results/logs/06_add_cluster_information_test.log

python -m cProfile -o ../results/profs/07_print_results_test.prof \
    07_print_results.py -db $DATABASE -o results_test \
    2>&1 | tee  ../results/logs/07_print_results_test.log