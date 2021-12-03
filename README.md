# bioproject

**1. Data acqusition**   
For this project we will be using the sequence of the *Ramazzottius varieornatus*, the YOKOZUNA-1 strain (sequenced in the University of Tokyo and named after the highest rank in professional sumo). The assembled genome here: http://kumamushi.org/data/YOKOZUNA-1.scaffolds.fa  
**1. Repeat masking**     
**2. Gene prediction**   
**3. Model training**   
  
These steps were skipped for now due to computational complexity.  

**4. Intro -  Functional annotation.**   
For this step downloaded provided augustus.whole.gff data.   
In order to extract proteins to fasta used perl script getAnnoFast.pl, which produced a file with proteins agustus.whole.fa. As found by record count there are 16435 proteins in a file.  
perl getAnnoFast.pl augustus.whole.gff
cat augustus.whole.fa | grep '>' | wc -l
**5. Physical localization**
/Users/aleksandrkovalenko/anaconda3/pkgs/blast-2.6.0-boost1.64_2/bin/makeblastdb -in augustus.whole.fasta -dbtype prot -out proteins_db

# Output: Adding sequences from FASTA; added 16435 sequences in 0.375326 seconds.
/Users/aleksandrkovalenko/anaconda3/pkgs/blast-2.6.0-boost1.64_2/bin/blastp -query peptides.fa -db proteins_db -out proteins.blastp -outfmt "6 qseqid sseqid evalue qcovs pident" -evalue 0.05 -task blastp-short

cat proteins.blastp | cut -f 2 > names.txt
# их как-то подозрительно мало, 44 белка

xargs samtools faidx augustus.whole.fasta < names.txt > localized_proteins.fasta
**6. Localization prediction**
There are still too many proteins to verify each of them experimentally. Now, we will try to predict where these proteins are found in the cell based on their sequences.  
- WoLF PSORT  (https://wolfpsort.hgc.jp/) predicts the subcellular localization of proteins.
g5641.t1 details extr: 31, lyso: 1
g12562.t1 details extr: 30, lyso: 2
g5641.t1 details extr: 31, lyso: 1
g5641.t1 details extr: 31, lyso: 1
g12562.t1 details extr: 30, lyso: 2
g5616.t1 details extr: 31, mito: 1
g13530.t1 details extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5
g13530.t1 details extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5
g14472.t1 details nucl: 28, plas: 2, cyto: 1, cysk: 1
g4106.t1 details E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1
g10513.t1 details nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1
g5616.t1 details extr: 31, mito: 1
g5616.t1 details extr: 31, mito: 1
g5502.t1 details extr: 31, lyso: 1
g5467.t1 details extr: 27, plas: 4, mito: 1
g15153.t1 details extr: 32
g5641.t1 details extr: 31, lyso: 1
g12562.t1 details extr: 30, lyso: 2
g5616.t1 details extr: 31, mito: 1
g5503.t1 details extr: 29, plas: 1, mito: 1, lyso: 1
g5502.t1 details extr: 31, lyso: 1
g5510.t1 details plas: 23, mito: 7, E.R.: 1, golg: 1
g12510.t1 details plas: 29, cyto: 3
g5237.t1 details plas: 24, mito: 8
g12510.t1 details plas: 29, cyto: 3
g12510.t1 details plas: 29, cyto: 3
g5641.t1 details extr: 31, lyso: 1
g12562.t1 details extr: 30, lyso: 2
g5641.t1 details extr: 31, lyso: 1
g15153.t1 details extr: 32
g5641.t1 details extr: 31, lyso: 1
g12562.t1 details extr: 30, lyso: 2
g5616.t1 details extr: 31, mito: 1
g5503.t1 details extr: 29, plas: 1, mito: 1, lyso: 1
g5502.t1 details extr: 31, lyso: 1
g4106.t1 details E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1
g5616.t1 details extr: 31, mito: 1
g5616.t1 details extr: 31, mito: 1
g5502.t1 details extr: 31, lyso: 1
g5467.t1 details extr: 27, plas: 4, mito: 1
g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
g12562.t1 details extr: 30, lyso: 2
g5237.t1 details plas: 24, mito: 8
# only proteins at least with some nuclear association (6 proteins are collected)  
cat WolFPSORT.txt | grep 'nucl' > WolFPSORT_nucl.txt  
g13530.t1 details extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5
g13530.t1 details extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5
g14472.t1 details nucl: 28, plas: 2, cyto: 1, cysk: 1
g10513.t1 details nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1
g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
- TargetP - 2.0. Also predicts localizatonn using terminal sequences. 

# TargetP-2.0	Organism: Non-Plant	Timestamp: 20211129133620
# ID	Prediction	OTHER	SP	mTP	CS Position
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g13530.t1	SP	0.116007	0.883840	0.000153	CS pos: 19-20. TIP-FT. Pr: 0.3552
g13530.t1	SP	0.116007	0.883840	0.000153	CS pos: 19-20. TIP-FT. Pr: 0.3552
g14472.t1	OTHER	0.999999	0.000001	0.000000	
g4106.t1	OTHER	0.729658	0.266917	0.003425	
g10513.t1	OTHER	0.999999	0.000001	0.000000	
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5502.t1	SP	0.001134	0.998823	0.000043	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5467.t1	SP	0.000096	0.999845	0.000059	CS pos: 16-17. ASA-GS. Pr: 0.6543
g15153.t1	SP	0.000014	0.999986	0.000000	CS pos: 16-17. AYA-AN. Pr: 0.8378
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5503.t1	SP	0.001222	0.998720	0.000058	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5502.t1	SP	0.001134	0.998823	0.000043	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5510.t1	OTHER	0.999108	0.000016	0.000876	
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g5237.t1	OTHER	0.999545	0.000345	0.000111	
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g15153.t1	SP	0.000014	0.999986	0.000000	CS pos: 16-17. AYA-AN. Pr: 0.8378
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5503.t1	SP	0.001222	0.998720	0.000058	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5502.t1	SP	0.001134	0.998823	0.000043	CS pos: 16-17. ASA-GS. Pr: 0.6833
g4106.t1	OTHER	0.729658	0.266917	0.003425	
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5502.t1	SP	0.001134	0.998823	0.000043	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5467.t1	SP	0.000096	0.999845	0.000059	CS pos: 16-17. ASA-GS. Pr: 0.6543
g15484.t1	OTHER	0.999980	0.000010	0.000010	
g15484.t1	OTHER	0.999980	0.000010	0.000010	
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g5237.t1	OTHER	0.999545	0.000345	0.000111	
# we are interested in 'OTHER peptides'

# ID	Prediction	OTHER	SP	mTP	CS Position
g14472.t1	OTHER	0.999999	0.000001	0.000000	
g4106.t1	OTHER	0.729658	0.266917	0.003425	
g10513.t1	OTHER	0.999999	0.000001	0.000000	
g5510.t1	OTHER	0.999108	0.000016	0.000876	
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g5237.t1	OTHER	0.999545	0.000345	0.000111	
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g4106.t1	OTHER	0.729658	0.266917	0.003425	
g15484.t1	OTHER	0.999980	0.000010	0.000010	
g15484.t1	OTHER	0.999980	0.000010	0.000010	
g5237.t1	OTHER	0.999545	0.000345	0.000111	
The overlap with WOLF:
- g14472.t1	OTHER	0.999999	0.000001	0.000000
- g10513.t1	OTHER	0.999999	0.000001	0.000000
- g15484.t1	OTHER	0.999980	0.000010	0.000010
**7. BLAST search**
Let's blast overlapping peptides:  
samtools faidx augustus.whole.fasta g14472.t1 g10513.t1 g15484.t1 > proteins_for_blast.fasta
- **g14472.t1**. no similarity found
- **g10513.t1**.
![image.png](attachment:image.png)
- **g15484.t1**. has 100 close matches. https://blast.ncbi.nlm.nih.gov/Blast.cgi It is though something like vacuolar protein sorting-associated protein

Let's blas some additional sequences: 
- **g13530.t1**. no similarity
**8. Pfam prediction**  
We will use HMMER (web-version, https://www.ebi.ac.uk/Tools/hmmer/) to search our protein sequences against a collection of profile-HMMs for different protein domains and motifs.  
Especially we will use hmmscan in Pfam database  
- **g14472.t1**. no hits
- **g10513.t1**. no hits
- **g15484.t1** https://www.ebi.ac.uk/Tools/hmmer/results/AE108FAC-5116-11EC-AB0D-7248F75AEC3D/score    
![image.png](attachment:image.png)
- **g13530.t1**. no hits
