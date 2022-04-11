#!/bin/bash


DATABASE="photorhabdus_xenorhabdus.db"
REF_GENOMES_DIR="actinobacteria_internal streptomyces_internal"
QUERY_GENOMES_DIR="photorhabdus xenorhabdus"


### Fist part on local machine, as the query genomes are not on the shared machine. Could be changed, if O drive can be connected to shared machine.
### to use the O drive locally (ubuntu subsystem on win) the drive needs to be mounted first: sudo mount -t drvfs O: /mnt/o

python -m cProfile -o ../results/profs/01_create_database.prof \
    01_create_database.py -db $DATABASE \
    2>&1 | tee ../results/logs/01_create_database.log

python -m cProfile -o ../results/profs/02_import_genomes_ref.prof \
    02_import_genomes.py -db $DATABASE -r $REF_GENOMES_DIR \
    2>&1 | tee ../results/logs/02_import_genomes_ref.log


### Second part on shared machine

python -m cProfile -o ../results/profs/02_import_genomes_query.prof \
    02_import_genomes.py -db $DATABASE -q $QUERY_GENOMES_DIR \
    2>&1 | tee ../results/logs/02_import_genomes_query.log

python -m cProfile -o ../results/profs/03_add_core_genome_indicator.prof \
    03_add_core_genome_indicator.py -db $DATABASE \
    2>&1 | tee ../results/logs/03_add_core_genome_indicator.log

python -m cProfile -o ../results/profs/04_find_approximate_pfamsequences.prof \
    04_find_approximate_pfamsequences.py -db $DATABASE -t 8 \
    2>&1 | tee ../results/logs/04_find_approximate_pfamsequences.log

python -m cProfile -o ../results/profs/05_cluster_hits.prof \
    05_cluster_hits.py -db $DATABASE -t 8 \
    2>&1 | tee ../results/logs/05_cluster_hits.log

python -m cProfile -o ../results/profs/06_add_cluster_information.prof \
    06_add_cluster_information.py -db $DATABASE \
    2>&1 | tee ../results/logs/06_add_cluster_information.log

python -m cProfile -o ../results/profs/06_print_results.prof \
    06_print_results.py -db $DATABASE \
    2>&1 | tee ../results/logs/07_print_results.log