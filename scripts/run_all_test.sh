#!/bin/bash

# to use the O drive locally (ubuntu subsystem on win) the drive needs to be mounted: sudo mount -t drvfs O: /mnt/o

DATABASE="database_test.db"

python -m cProfile -o profs/01_create_database_test.prof 01_create_database.py -db $DATABASE > logs/01_create_database_test.log 2>&1
python -m cProfile -o profs/02_import_genomes_test.prof 02_import_genomes.py -db $DATABASE --test > logs/02_import_genomes_test.log 2>&1
python -m cProfile -o profs/03_add_core_genome_indicator_test.prof 03_add_core_genome_indicator.py -db $DATABASE > logs/03_add_core_genome_indicator_test.log 2>&1
python -m cProfile -o profs/04_find_approximate_pfamsequences_test.prof 04_find_approximate_pfamsequences.py -db $DATABASE -t 2 > logs/04_find_approximate_pfamsequences_test.log 2>&1
python -m cProfile -o profs/05_cluster_hits_test.prof 05_cluster_hits.py -db $DATABASE -t 2 > logs/05_cluster_hits_test.log 2>&1
python -m cProfile -o profs/06_add_cluster_information_test.prof 06_add_cluster_information.py -db $DATABASE > logs/06_add_cluster_information_test.log 2>&1
python -m cProfile -o profs/07_print_results_test.prof 06_print_results.py -db $DATABASE -o results_test > logs/07_print_results_test.log 2>&1