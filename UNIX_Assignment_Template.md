
#UNIX Assignment

##Data Inspection

###Attributes of `fang_et_al_genotypes`

```
# to know the  size
ls -lh fang_et_al_genotypes.txt

# to know the file type
file fang_et_al_genotypes.txt

# count number of lines 
wc -l fang_et_al_genotypes.txt 

# count number of columns in the header
head -n 1 fang_et_al_genotypes.txt | awk -F'\t' '{print NF}'

# inspect the first few columns 
cut -f1-10 fang_et_al_genotypes.txt | head -n 5

# check the groups and the count of each
cut -f3 fang_et_al_genotypes.txt | tail -n +2 | sort | uniq -c | sort -nr | head

# check if the Sample_IDs is unique 
cut -f1 fang_et_al_genotypes.txt | tail -n +2 | sort | uniq -d | head
```

By inspecting this file I learned that:

1. The file size is 11M 
2. it is ASCII text with very long lines.
3. it has 2783 lines (including 1 header line)
4. It has 986 columns in header line (tab delimited)
5. The first 3 columns are metadata (Sample_ID, JG_OTU, Group).
6. The rest of the columns are SNP genotype columns such as abph1.20, abph1.22, ae1.3
7. There are 10 groups. Here is a list of them along with the count of each.
1256 ZMMLR
 900 ZMPBA
 290 ZMMIL
  75 ZMXCH
  69 ZMXCP
  41 ZMPIL
  34 ZMPJA
  27 ZMMMR
  22 TRIPS
  17 ZLUXR
8. No duplication in Sample_ID column, so they are all unique.





###Attributes of `snp_position.txt`

```
# to know the  size
ls -lh snp_position.txt

# to know the file type
file snp_position.txt

# count number of lines 
wc -l  snp_position.txt

# count number of columns in the header
head -n 1 snp_position.txt         | awk -F'\t' '{print NF}'

# inspect the first few columns
cut -f1-6 snp_position.txt | head -n 10

# check the chromosomes  and the count of each
cut -f3 snp_position.txt | tail -n +2 | sort | uniq -c | sort -nr
```

By inspecting this file I learned that:

1.  The file size is 81K
2.it is ASCII text
3.  It has 984 lines (including the header)
4. It has 15 columns in header line (tab delimited)
5. Here is a list of the chromosomes (10 chromosomes) included and their occurences

 155 1
 127 2
 122 5
 107 3
  97 7
  91 4
  76 6
  62 8
  60 9
  53 10
Which means Chromosome 1 has 155 SNPs, chr2 127, chr5 122, chr3 107, chr7 97, chr4 91, chr6 76, chr8 62, chr9 60, chr10 53.
6. There are 27 SNPs with the label unknown and 6 SNPs with the label multiple


##Data Processing

###Maize Data

```
awk -F'\t' 'NR==1 || $3=="ZMMIL" || $3=="ZMMLR" || $3=="ZMMMR"' fang_et_al_genotypes.txt > maize_genotypes.txt
wc -l maize_genotypes.txt 
awk -f transpose.awk maize_genotypes.txt > maize_transposed.txt
head -n 6 maize_transposed.txt
awk -F'\t' 'NR==1{c=NF} END{print NR " rows x " c " cols"}' maize_genotypes.txt
awk -F'\t' 'NR==1{c=NF} END{print NR " rows x " c " cols"}' maize_transposed.txt
head -n 1 maize_transposed.txt | sed 's/^Sample_ID/SNP_ID\tChromosome\tPosition/' > maize_header.txt
awk 'NR>3' maize_transposed.txt | sort -k1,1 > maize_snps_sorted.txt
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1,$3,$4}' snp_position.txt | sort -k1,1 > positions_sorted.txt
join -1 1 -2 1 -t $'\t' positions_sorted.txt maize_snps_sorted.txt > maize_joined.txt
awk -F'\t' 'BEGIN{OFS="\t"} $2=="unknown" || $3=="unknown" {print}' maize_joined.txt > maize_unknown.txt
awk -F'\t' 'BEGIN{OFS="\t"} $2=="multiple" {print}' maize_joined.txt > maize_multiple.txt
awk -F'\t' 'BEGIN{OFS="\t"} $2!="unknown" && $2!="multiple" && $3!="unknown" && $3!="" {print}' maize_joined.txt > maize_single.txt
wc -l maize_joined.txt maize_unknown.txt maize_multiple.txt maize_single.txt

sed 's#\?/\?#?#g' maize_single.txt > maize_single_q.txt
sed 's#\?/\?#-#g' maize_single.txt > maize_single_dash.txt
for chr in {1..10}; do
  awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"} $2==c {print}' maize_single_q.txt | sort -t $'\t' -k3,3n > maize_chr${chr}_pos_asc.txt
  awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"} $2==c {print}' maize_single_dash.txt  | sort -t $'\t' -k3,3nr > maize_chr${chr}_pos_desc.txt
done



```

1- Separated the maize genotype  (ZMMIL, ZMMLR, ZMMMR)
2- check if the maize split is correct. It is correct because the code output was 1574 for and we know that ZMMLR = 1256, ZMMIL = 290, and ZMMMR = 27, therefore, the Sum = 1256 + 290 + 27 = 1573 data lines (+1 line header it becomes 1574)
3- Use the given transpose script to transpose the maize genotype file 
4- inspecting the first few lines in the transposed file. It looks alright because the rows now start with the column header in the original file (i.e. Sample_ID, JG_OTU, Group).
5- Checking if the transpose was correct using the dimensions. It was correct because the transposed maize file was 986 rows x 1574 cols cols and the original maize file has 1574 rows x 986 cols
6- Since all our files should be formatted such that the first column is "SNP_ID", the second column is "Chromosome", the third column is "Position", and subsequent columns are genotype data, I used the given code to create the required header for final outputs. Itconverts: Sample_ID ZDP_... into SNP_ID Chromosome Position ZDP_... (i.e. adds the Chromosome as column 2 and Position as column 3 and the rest of the column change their order accordingly.
7- We needed to remove the first 3 rows which contain metadata before we were able to join 
8- Rearranged the positions file so that SNP_ID stays in column 1, Chromosome becomes in column 2 and Position becomes in column 3 then sorted according to SNP_ID.
9- we joined the positions with maize genotypes file
10- Split the  UNKNOWN and MULTIPLE for maize
11-check the numbers of duplicates and unknowns to see if they agree with what we inspected from the beginning(e.g. 6 duplicates, 27 unknowns, etc,).  They agree.
12- Replaced missing genotype encoding using the symbol "?" for ascending order and the symbol "-" for descending order.
13- We looped chromosomes 1–10, then sorted them by position to create 10 chromosome files ascending + 10 descending for maize


###Teosinte Data

```
awk -F'\t' 'NR==1 || $3=="ZMPBA" || $3=="ZMPIL" || $3=="ZMPJA"' fang_et_al_genotypes.txt > teosinte_genotypes.txt
wc -l teosinte_genotypes.txt
awk -f transpose.awk teosinte_genotypes.txt > teosinte_transposed.txt
head -n 6 teosinte_transposed.txt
awk -F'\t' 'NR==1{c=NF} END{print NR " rows x " c " cols"}' teosinte_genotypes.txt
awk -F'\t' 'NR==1{c=NF} END{print NR " rows x " c " cols"}' teosinte_transposed.txt  
head -n 1 teosinte_transposed.txt | sed 's/^Sample_ID/SNP_ID\tChromosome\tPosition/' > teosinte_header.txt
awk 'NR>3' teosinte_transposed.txt | sort -k1,1 > teosinte_snps_sorted.txt
join -1 1 -2 1 -t $'\t' positions_sorted.txt teosinte_snps_sorted.txt > teosinte_joined.txt
awk -F'\t' 'BEGIN{OFS="\t"} $2=="unknown" || $3=="unknown" {print}' teosinte_joined.txt > teosinte_unknown.txt
awk -F'\t' 'BEGIN{OFS="\t"} $2=="multiple" {print}' teosinte_joined.txt > teosinte_multiple.txt
awk -F'\t' 'BEGIN{OFS="\t"} $2!="unknown" && $2!="multiple" && $3!="unknown" && $3!="" {print}' teosinte_joined.txt > teosinte_single.txt
wc -l teosinte_joined.txt teosinte_unknown.txt teosinte_multiple.txt teosinte_single.txt

sed 's#\?/\?#?#g' teosinte_single.txt > teosinte_single_q.txt
sed 's#\?/\?#-#g' teosinte_single.txt > teosinte_single_dash.txt
for chr in {1..10}; do
  awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"} $2==c {print}' teosinte_single_q.txt | sort -t $'\t' -k3,3n > teosinte_chr${chr}_pos_asc.txt
  awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"} $2==c {print}' teosinte_single_dash.txt  | sort -t $'\t' -k3,3nr > teosinte_chr${chr}_pos_desc.txt
done


```

1- Separate the teosinte genotype (ZMPBA, ZMPIL, ZMPJA)
2- check if the teosinte split is correct. It is correct because the code output was 976 for and we know that ZMPBA = 900, ZMPIL = 41, and ZMPJA = 34, therefore, the Sum = 900 + 41 + 34 = 975 data lines (+1 line header it becomes 976)
3- Use the given transpose script to transpose the teosinte genotype file 
4- inspecting the first few lines in the transposed file. It looks alright because the rows now start with the column header in the original file (i.e. Sample_ID, JG_OTU, Group).
5- Checking if the transpose was correct using the dimensions. It was correct because the transposed teosinte file was 986 rows x 976 cols and the original teosinte file has 976 rows x 986 cols.
6- Since all our files should be formatted such that the first column is "SNP_ID", the second column is "Chromosome", the third column is "Position", and subsequent columns are genotype data, I used the given code to create the required header for final outputs. Itconverts: Sample_ID ZDP_... into SNP_ID Chromosome Position ZDP_... (i.e. adds the Chromosome as column 2 and Position as column 3 and the rest of the column change their order accordingly.

7- We needed to remove the first 3 rows which contain metadata before we were able to join 
8-  we joined the positions with teosinte genotypes file
9- Split the  UNKNOWN and MULTIPLE for teosinte
10- checked the numbers of duplicates and unknowns to see if they agree with what we inspected from the beginning(e.g. 6 duplicates, 27 unknowns, etc,).  They agree.
11- Replaced missing genotype encoding using the symbol "?" for ascending order and the symbol "-" for descending order.
12- We looped chromosomes 1–10, then sorted them by position to create 10 chromosome files ascending + 10 descending for teosinte




