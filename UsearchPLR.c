pwd #print working directory#
cd ~/Documents/Pipelines/Community_Usearch  #Change working directory to Docs then pipelines then the Usearch pipeline where Usearch is located#
chmod a+rx UsearchExec # run Usearch executive#

#--Common UMGC primer sets for reference--#
#V's are areas on ITS. Different V's have different accuracies
#Primer set for V1-V3
AGAGTTTGATCMTGGCTCAG
ATTACCGCGGCTGCTGTGG

#primer seqs for V4 region
515F 5'-GTGCCAGCMGCCGCGGTAA
806R 5'-GGACTACHVHHHTWTCTAAT

#primer seqs for V3-V4 region
V3F:5'-CCTACGGGAGGCAGCAG
806R:5'-GGACTACHVGGGTWTCTAAT

#Primer set for ITS1
CTTGGTCATTTAGAGGAAGTAA
GCTGCGTTCTTCATCGATGC

#Primer set for ITS2
TCGATGAAGAACGCAGCG
TCCTCCGCTTATTGATATGC

###Setting Trimming Parameters###
./usearch11.exe -fastx_truncate 1051_S132_L001_R1_001.fastq -stripleft 22 -stripright 0 -fastqout 1051_S132_L001_R1_001.fastq.trim ## Stripping setup for the forward read 1051 is the first sample R1 is forward read

./usearch11.exe -fastx_truncate 1051_S132_L001_R2_001.fastq -stripleft 20 -stripright 20 -fastqout 1051_S132_L001_R2_001.fastq.trim ## Stripping setup for the reverse read (R2)##

./usearch11.exe -fastq_mergepairs 1051_S132_L001_R1_001.fastq.trim -relabel '@' -fastq_maxdiffs 5  -fastqout 1051_S132_L001_R1_001.fastq.trim.merged.fq  ##This merges the two reads but I am not sure why the R2 file isn't used ASK DAN

#Once desired trimming paramters are defined, use bash loops to apply commands to each file.

#loop for trimming low-quality ends and read pairing
for F in *R1*.fastq; do
	./usearch11.exe -fastx_truncate $F -stripleft 22 -stripright 0 -fastqout $F.trim.fq
done
####HOW COOL YOU NOW KNOW THE SYNTAX OF A FORLOOP IN SHELLSCRIPT!!! for F sets the variable, *R1*.fastq may mean includes R1 in fastq file. so then in the body of the forloop command $F replaces file names. 
for F in *R2*.fastq; do
./usearch11.exe -fastx_truncate $F -stripleft 20 -stripright 20 -fastqout $F.trim.fq
done

#loop for merging reads
for F in *R1_001*trim.fq; do
	./usearch11.exe -fastq_mergepairs $F -relabel @ -fastq_maxdiffs 5 -fastqout $F.merged.fq
done
####So there is something that -fastq_mergepairs does that I'm not super sure about but I think that it will read R1 and R2 and match them. I'll look at the help section and update this note!!!!

#Convert merged reads to fasta format and combine into a single pool of reads. 
#These reads will be mapped to denoised sequences later on to generate an OTU table.
for f in *merged.fq; do
	./usearch11.exe -fastq_filter $f -fastaout $f.fa
done

cat *merged.fq.fa > reads.fa
rm *merged.fq.fa #removes intermediate merged fasta files

##################
#Step 2: Generate high-quality subset of reads for denoising to generate amplicon sequence variants (aka Zotus, OTUs, ESVs, etc)
#################


#Quality filter with maximum expected error to generate high-quality read set for denoising. It may also
#be useful to change other paramenters (eg. minimum or maximum length) depending on the expected amplicon size, etc 
for F in *merged.fq; do
./usearch11.exe -fastq_filter $F -fastq_maxee 1 -fastq_minlen 200 -fastaout $F.filtered.fa
done
##### Ok so we filter out everything with an error of 1 or more which I think means the  reverse didn't match at one or more places "-fastq_maxee 1" 
#####and the length of the read hast to be at least 200bp "-fastq_minlen 200"
##### output as a fasta, title refers to the process of removing more than 1bp error rate "filtered" in "filtered.fa"

#combine filtered reads
cat *filtered.fa > filtered.all.fa #DAN!!!!
rm *filtered.fa #removes intermediate filtered fasta files
#####I'm not sure if cat means concatenate or not I will look at the literature Usearch provides

#dereplicate sequences (reduce to unique sequences only). Vsearch is sometimes necessary here if usearch memory caps are reached.
./usearch11.exe -fastx_uniques filtered.all.fa -fastaout uniques.fa -sizeout -relabel Uniq
#### The language here is self explanatory. I am not sure what it deems unique since maybe there is some cross over in reads I will ask DAN!!!!!!!

#./vsearch.exe -derep_fulllength filtered.all.fa -output uniques.fa -sizeout -relabel Uniq

#denoise sequences with unoise3 algorithm. Vsearch is sometimes necessary here if usearch memory caps are reached.
./usearch11.exe -unoise3 uniques.fa -zotus zotus.fa
####I forgot what the denoise step actually does LOOK AT NOTES!!!!! IF CANT FIND ASK DAN!!!!!

#denoise with vsearch, if necessary
#./vsearch.exe -cluster_unoise uniques.fa consout - zotus.fa
#additional chimera detection with vsearch -- necessary if using denoising is performed with vsearch
#./vsearch.exe -uchime3_denovo zotus.fa --nonchimeras zotus_nochim.fa --chimeras zotus_chimeras.fa


####################################################
#Step 3: Map merged read pool to representative ASVs. 
####################################################

#create an OTU table by mapping reads back to zOTUs. Vsearch is sometimes necessary here.
./usearch11.exe -otutab reads.fa -zotus zotus.fa -otutabout zotutab.txt -sample_delim .
####Looks like there's a lot going on here. Read the Usearch lit to understand more

#./vsearch.exe -usearch_global reads.fa --db zotus_nochim.fa --id 0.97 --otutabout zotutab.txt


###################################################
#Step 4: Assign taxonomy to ASVs. 
 #Note this is one of many taxonomy assignment methods, and other approaches/databases may produce different results
####################################################

#classify zotus using SINTAX k-mer based approach
./usearch11.exe -sintax zotus.fa -db utax_reference_dataset_10.05.2021.fasta -tabbedout zotus.sintax -strand both -sintax_cutoff 0.8 

#####BAM here is the stept that matches OTU to your data set. Should read the literature on this as well. I have a feeling RNAseq prep will be a similar process but we'll have to look at it down the line.

#classify zotus using naieve bayesian classifier
./usearch11.exe -nbc_tax zotus.fa -db rdp_16s_v18.fa -tabbedout zotus.nbc -strand plus


##The important output files from this pipeline are:
# 1. 'zotus.fa' a fasta of ASV sequences
# 2. 'zotutab.txt' an OTU table of sequence counts/sample
# 3. 'zotus.sintax' the taxonomy assignment of ASVs 




