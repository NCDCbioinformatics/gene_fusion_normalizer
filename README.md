# gene_fusion_normalizer
<img width="2554" height="915" alt="image" src="https://github.com/user-attachments/assets/75577907-821e-4f22-b301-fd749f971ac4" />

#install 

unzip gene_fusion_normalizer_0.2.1 \
cd gene_fusion_normalizer \
pip install -e . \



#run --col Automatically recognizes fusion col when no input is provided

  gene-fusion-normalizer "/path/gene_fusion_normalizer/gene_split_test.xlsx" \
  --col "fusion" \
  --gtf "/path/Homo_sapiens.GRCh38.110.gtf.gz" \
  --hgnc "/path/hgnc_complete_set.txt" \
  --explode \
  -o "/out_path/fusion_mapped.xlsx"

  
# Feature Summary
- Standardize and separate fusion gene names (geneA-geneB) written in various ways
- Basically, fusion genes are divided based on "-", but if there are multiple "-" \
  => Split all cases and check whether the corresponding gene name is known \
  => If it is an intergenic or RNA name, enter it as is
 

