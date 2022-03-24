#!/bin/bash


DATABASE="database.db"

### First 2 scripts were run locally, because I dont know how to connect the O drive to the server. To use the O drive locally (ubuntu subsystem on win) the drive needs to be mounted: sudo mount -t drvfs O: /mnt/o

python -m cProfile -o ../results/profs/01_create_database.prof \
    01_create_database.py -db $DATABASE \
    > ../results/logs/01_create_database.log 2>&1

python -m cProfile -o ../results/profs/02_import_genomes.prof \
    02_import_genomes.py -db $DATABASE \
    > ../results/logs/02_import_genomes.log 2>&1

### Next scripts were run on the server

python -m cProfile -o ../results/profs/03_add_core_genome_indicator.prof \
    03_add_core_genome_indicator.py -db $DATABASE \
    > ../results/logs/03_add_core_genome_indicator.log 2>&1

python -m cProfile -o ../results/profs/04_find_approximate_pfamsequences.prof \
    04_find_approximate_pfamsequences.py -db $DATABASE -t 12 \
    > ../results/logs/04_find_approximate_pfamsequences.log 2>&1

python -m cProfile -o ../results/profs/05_cluster_hits.prof \
    05_cluster_hits.py -db $DATABASE -t 8 \
    > ../results/logs/05_cluster_hits.log 2>&1

python -m cProfile -o ../results/profs/06_add_cluster_information.prof \
    06_add_cluster_information.py -db $DATABASE \
    > ../results/logs/06_add_cluster_information.log 2>&1

python -m cProfile -o ../results/profs/06_print_results.prof \
    06_print_results.py -db $DATABASE \
    > ../results/logs/07_print_results.log 2>&1