#!/bin/bash

# USAGE: bash magma.sh <working dir> <binary file prefix> <gene location file> <path to sample dir> <phenotype file> <gene set file> <output prefix>

wdir=$1
bfile=$2
geneloc=$3
sample_dir=$4
phenofile=$5
genesetfile=$6
output=$7

cd ${wdir}

#dirlist=(subsets annotates magma_gene magma_geneset)
#for dir in  ${dirlist[@]}
#do
#    if [[ -d ${dir} ]]
#    then
#        echo "Folder ${dir} exists, making new folder..."
#        rm -r ${dir}
#        mkdir ${dir}
#    else
#        echo Making ${dir}
#        mkdir ${dir}
#    fi
#done
#
#echo Annotating SNPs to genes...  | sed G
#
#./magma --annotate window=2,0.5 --snp-loc ${bfile}.bim --gene-loc ${geneloc} --out ./annotates/annotates
#
#echo Starting geneset analysis | sed G
#
#for i in {1..100};
#do
#    echo ...Iteration 1
#	###Removing training samples
#	./plink --silent --bfile ${bfile} --keep ${sample_dir}/samples_${i}.txt --make-bed --out ./subsets/${output}_${i}
#
#	echo ...Running gene-level analysis
#    # Gene-level analysis for 1 bootstrapped subset
#	./magma --bfile ./subsets/${output}_${i} --gene-annot ./annotates/annotates.genes.annot --pheno file=${phenofile} --out ./magma_gene/${output}_${i}
#  
#	echo ...Running gene set-level analysis | sed G
#    # Gene-set level analysis for 1 bootstrapped subset
#	./magma --gene-results ./magma_gene/${output}_${i}.genes.raw --model alpha=0.3 --settings gene-info --set-annot ${genesetfile} --out ./magma_geneset/${output}_${i}
#
#done
#
echo Extract SNPs associated to significant genes in gene-level analysis | sed G
# Extract and aggregate significant genes 
Rscript extract_sig_genes.R $wdir $output

if [ -f ${output}_associated_SNPs.txt ]; then
    nr=$(wc -l < ${output}_associated_SNPs.txt)
    if [ "$nr" -eq 1 ]; then
        echo "...There is only 1 associated SNP" | sed G
        cp ${output}_associated_SNPs.txt ${output}_associated_LDprunned_SNPs.txt
    else 
        echo "...Get LDprunned SNPs" | sed G
        # Extract reduced bed file for prunning 
        ./plink --silent -bfile ${bfile} --extract ${output}_associated_SNPs.txt --keep ${sample_dir}/samples_fs.txt --make-bed --out tmp
        # Look for LD SNPs
        ./plink --bfile tmp --indep 50 5 2 
        cp plink.prune.in ${output}_associated_LDprunned_SNPs.txt
    fi
else
    echo There is no significant genes and SNPs to extract | sed G
fi

echo Extract leading-edge SNPs associated with significant pathways in geneset-level analysis | sed G
# Extract entire feature selection set
./plink --silent --bfile ${bfile} --keep ${sample_dir}/samples_fs.txt --make-bed --out fs_${output}

# Gene analysis for entire feature selection set
./magma --bfile ./fs_${output} --gene-annot ./annotates/annotates.genes.annot --pheno file=${phenofile} --out ./magma_gene/${output}

# Final GSEA to extract the leading edge SNPs
Rscript final_gsea.R ${output} $wdir

if [ -f ${output}_leadingEdge_SNPs.txt ]; then
    nr=$(wc -l < ${output}_leadingEdge_SNPs.txt)
    if [ "$nr" -eq 1 ]; then
        echo "...There is only 1 leading edge SNP" | sed G
        cp ${output}_leadingEdge_SNPs.txt ${output}_leadingEdge_LDprunned_SNPs.txt
    else 
        echo "...Get LDprunned SNPs" | sed G
        # Extract redcued bed file for prunning 
        ./plink --silent -bfile ${bfile} --extract ${output}_leadingEdge_SNPs.txt --keep ${sample_dir}/samples_fs.txt --make-bed --out tmp

        # Look for LD SNPs
        ./plink --bfile tmp --indep 50 5 2 
        cp plink.prune.in ${output}_leadingEdge_LDprunned_SNPs.txt
    fi
else
    echo There is no significant genes and SNPs to extract | sed G
fi


echo Extract final binary file set | sed G
cat ${output}_associated_LDprunned_SNPs.txt ${output}_leadingEdge_LDprunned_SNPs.txt | sort | uniq > final_SNPs.txt
# Make final bed file
./plink --silent -bfile ${bfile} --extract final_SNPs.txt --make-bed --out ${output}_genomics_selected

