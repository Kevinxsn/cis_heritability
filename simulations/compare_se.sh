for i in {1..100}; do
    python simulation_small0.2_10.py --geno_prefix /projects/ps-renlab2/sux002/data/test_gene_region >> results.txt
done

for i in {1..100}; do
    python simulation_small0.4_10.py --geno_prefix /projects/ps-renlab2/sux002/data/test_gene_region >> results.txt
done