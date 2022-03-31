#!/bin/bash

# to use the O drive locally (ubuntu subsystem on win) the drive needs to be mounted first: sudo mount -t drvfs O: /mnt/o

DATABASE="database_test.db"

python -m cProfile -o ../results/profs/01_create_database_test.prof \
    01_create_database.py -db $DATABASE \
    2>&1 | tee ../results/logs/01_create_database_test.log 

python -m cProfile -o ../results/profs/02_import_genomes_test.prof \
    02_import_genomes.py -db $DATABASE --test \
    2>&1 | tee  ../results/logs/02_import_genomes_test.log

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