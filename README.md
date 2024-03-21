![](static/protist.png)
![](static/fungi.tiff)
# The ProFun (Protist-Fungi) DB!
## By: Samantha Gleich & Syrena Whitner
## Updated: March 21, 2024
The Marine Microbial Eukaryote Transcriptome Sequencing Project (MMETSP) database is often used to annotate environmental sequence data; however, this databse does not contain references for a key group of marine protists, the MArine STramenopiles (MAST). Additionally, there is only one fungal group included in the MMETSP in its current form. Here, we show how to incorporate MAST single amplified genomes (SAGs) and fungal protein sequences into a EUKulele database (https://eukulele.readthedocs.io/en/latest/). This database can then be used to annotate environmental sequence data. 
<br>
<br>
This repository will require the use of [TransDecoder](https://github.com/TransDecoder/TransDecoder) and [EUKulele](https://github.com/AlexanderLabWHOI/EUKulele).


## Download SAGs and fungal protein sequences.
SAG data were obtained from Labarre et al. (2021); https://doi.org/10.1038/s41396-020-00885-8
<br>
[Access the MAST SAG data here](https://figshare.com/articles/dataset/Co-assembly/12430790?backTo=/collections/Comparative_genomics_reveals_new_functional_insights_in_uncultured_MAST_species/5008046)
<br>
<br>
Fungal reference genomes can be obtained from the MycoCosm database. We pulled all of the mycocosm protein files using the get_jgi_genomes GitHub repo. 
<br>
[Access the MycoCosm database here](https://mycocosm.jgi.doe.gov/mycocosm/home)
<br>
[Access the get_jgi_genomes Github repo here](https://github.com/guyleonard/get_jgi_genomes)

## Run TransDecoder on MAST SAGs to obtain putative protein-coding sequences.
```
./TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t Assembly.fasta
```
## Add Source ID to all MAST protein file sequence headers.
The SourceID needs to be in the header line (i.e., line containig '>') of each contig in each assembly. The SourceID you include here will need to match a taxonomy file that we will create next. 
```
sed 's/>.*/& \/SOURCE_ID=NameOfOrganism/‘ MAST_1.fasta > new_MAST_1.fa
rm MAST_1.fasta
```
## Add Source ID to all fungi protein file sequence headers.
Because we are working with > 1,000 fungal protein files, we added the Source ID to the sequence headers using a custom python script. See script add_source.py in this repository. 
```
python ./add_source.py
```

## Concatenate all new (non-MMETSP) protein files together. Concatenate new protein files and the MMETSP protein file. 
Each .fasta file here corresponds to a different reference that is being added to the MMETSP. Make sure that the files being concatenated all have the SOURCE_ID added in the previous step(s). 
```
cat new_*.fasta > all_fungal_proteins.fa
cat new_MAST*.fasta > all_MAST_proteins.fa
cat all_fungal_proteins.fa all_MAST_proteins.fa mmetsp.fa > fungi_mast_mmetsp.fa
```
## Make taxonomy table with taxonomy information for all of the new assemblies that are being added to the MMETSP. 
Make a tab separated taxonomy table that is in the same format as the EUKulele taxonomy table. Example: 
<br>
<br>
&emsp; Unnamed: 0 &emsp; Domain &emsp; Supergroup &emsp; Division &emsp; Class &emsp; Order &emsp; Family &emsp; Genus &emsp; Species &emsp; Source_ID
<br>
0 &emsp; 0 &emsp; Eukaryota &emsp; Opisthokonta &emsp; Ascomycota &emsp; Sordariomycetes &emsp; Microascales &emsp; Halosphaeriac &emsp; Corollospora &emsp; Corollospora_maritima &emsp; Corollospora
<br>
1 &emsp; 1 &emsp; Eukaryota &emsp; Opisthokonta &emsp; Basidiomycota &emsp; Agaricomycetes &emsp; Agaricales &emsp; Niaceae &emsp; Digitatispora &emsp; Digitatispora_marina &emsp; Digitatispora
<br>
<br>
The Source_ID columns must match the SOURCE_ID labels that were added to the .fa files above. Save this file as a .txt file (e.g., fungi_tax.txt) and remove the header row (i.e., Unnamed: 0 &emsp; Domain &emsp;) from the file for the next step. 

## Concatenate new taxonomy table with MMETSP taxonomy table.
```
cat mmetsp_tax.txt fungi_tax.txt mast_tax.txt > ProFun_Taxonomy.txt
```
## Run create_protein_table.py script in EUKulele to make a new database using the concatenated protein file and the concatenated taxonomy file.
```
create_protein_table.py --infile_peptide ProFun_Pep.fa --infile_taxonomy fungi_mast_mmmetsp_tax.txt --outfile_json ProFun_Map.json --output ProFun_Taxonomy.txt --delim "/"
```
## Run EUKulele with custom Pro(tist) Fun(gi) database. The previous EUKulele command (create_protein_table.py) will produce 2 new database files (fungi_mast_mmetsp_db.txt and fungi_mast_mmetsp_db.json).
```
EUKulele -m mets --sample_dir /path/to/directory --out_dir /path/to/directory/eukulele_out --reference_dir /path/to/directory/ --ref_fasta ProFun_Pep.fa --n_ext cds --tax_table ProFun_Taxonomy.txt --protein_map ProFun_Map.json
```
This EUKulele command will annotate your metatranscriptome assembly contigs using the custom Pro Fun database!
