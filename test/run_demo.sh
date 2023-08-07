PYTHONPATH=$PYTHONPATH:..

python ../waafle/waafle_search.py \
       ../demo/input/demo_contigs.fna \
       ../demo/input/demo_waafledb/demo_waafledb \

python ../waafle/waafle_genecaller.py \
       demo_contigs.blastout

python ../waafle/waafle_orgscorer.py \
       ../demo/input/demo_contigs.fna \
       demo_contigs.blastout \
       demo_contigs.gff \
       ../demo/input/demo_taxonomy.tsv \
       --write-details \
       --basename test \ 
