# RVTT
This code is associated with the following manuscript. If you use any part of the source code, please cite us:

"Pathway-based Rare Variant Burden Analysis Identifies a Role for the Complement System in An Extreme Phenotype of Sepsis with Coagulopathy"

Pavan K. Bendapudi, Sumaiya Nazeen, Justine Ryu, Onuralp Söylemez, Betty Rouaisnel, Meaghan Colling, Bryce Pasko, Alissa Robbins, Michael Bouzinier, Lindsay Tomczak, Lindsay Collier, Sanjay Ram, Agnes Toth-Petroczy, Joel Krier, Elizabeth Fieg, Walter H. Dzik, James C. Hudspeth, Olga Pozdnyakova, Valentina Nardi, Richard Maas, Shamil Sunyaev, and Julie-Aurore Losman.

## Requirements:
  - Python >= 3.6

## Instructions:
  - The input vcf file must be quality controlled as described in the manuscript and then annotated with SnpEff, SnpSIFT, and dbNSFP database as follows:
    - Annotation with SnpEff and SnpSift:
      
      `java -Xmx40g -jar $SNPEFF_BIN <SNPEFF_DATABASE> <QCED_VCF_GZ> > <SNPEFF_ANNOTATED_VCF>`
      
      `java -Xmx40g -jar $SNPSIFT_BIN dbnsfp -v -db <SNPSIFT_DATABASE> -f hg19_chr,"hg19_pos(1-based)", aaref, aaalt, aapos, genename, Ancestral_allele, SIFT_pred, Polyphen2_HVAR_pred, CADD_phred, gnomAD_exomes_NFE_AF, gnomAD_exomes_POPMAX_AF, gnomAD_genomes_NFE_AF, gnomAD_genomes_POPMAX_AF, clinvar_hgvs, Interpro_domain, GTEx_V7_gene, GTEx_V7_tissue -a -m <SNPEFF_ANNOTATED_VCF>  > <SNPSIFT_ANNOTATED_VCF>`
      
      `bgzip -c <SNPSIFT_ANNOTATED_VCF> > <SNPSIFT_ANNOTATED_VCF_GZ>`
      
      `tabix <SNPSIFT_ANNOTATED_VCF_GZ>`
      
  - The annotated and filtered gzipped vcf file can be preprocessed using the preprocess_gzvcf.py program
    
    Command: `python3 preprocess_gzvcf.py filtered_input.vcf.gz famfile mincutoff maxcutoff outprefix`
    
  - RVTT with fixed threshold can be run as follows:
    
    Command: `python3 rvtt_fixed_threshold.py <output of preprocess_gzvcf.py> <genelist> <famfile> <cutoff> <N> <seed> <outfile>`
 
  - RVTT with fixed threshold can be run as follows:
    
    Command: `python3 rvtt_variable_threshold.py <output of preprocess_gzvcf.py> <genelist> <famfile> <cutoff> <N> <seed> <outfile>`
