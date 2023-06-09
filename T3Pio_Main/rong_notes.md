# Rong Jin's notes on running Sean's T3pio pipeline


### 1.  install required packages  

I put together the t3pio.yaml, so you can create a conda env with it. For that, you will run the following:   
  >   1. `conda env create -n t3pio -f t3pio.yaml` (or `mamba env create -n t3pio -f t3pio.yaml` instead for speed; CDC users have mamba available once they module load conda)   
  >   2.  `conda activate t3pio`  



### 2.  run the pipeline  

You can execute the moduleFile.py script by typing:  
`python3 moduleFile.py num_gbk_files perc_of_inclusion boulderFile /path-to-gbk-file/ output_folder num_cpu`  

For example:  
`python3 moduleFile.py 19 1 boulderFile /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/gbks/ 19gbks_output1 20`

### 3.  MUSCLE output  
It turns out in Sean's script, MUSCLE grabs all available CPU and sends a ton of output to screen. So I instead changed the script and set threads to 2 and suppress output to screen.
`return_code = subprocess.check_output(['muscle','-align',multifastaFile,'-output',muscleFile+'.muscle', '-threads', '2', '-quiet'])`  


### 4.  The MasterPrimerFile   

 **located at /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/MasterPrimerFile**  Sean's note:  
1. List of the identified shared primers from the blastANI approach and stoolbug (JulieBugs) testing
2. From the pipeline version 0.3

>  **It has total 7600 primers, but how did we create it ? Using which collection of genomes ?**  

>  **2nd question are our 2461 primers a subset of these 7600 primers**  
>  -- Yes, I confirmed this by running this command  
>  `grep -xFf primer_pairs_3columns_rearranged MasterPrimerFile | wc -l` in the `/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data` directory, and it's 2461, all those 2461 primers can be found in the MasterPrimerFile.    
>  `x: entire line matches a pattern`  
>  `F: literal string, no regex`  
>  `f: patterns in a file`  

### 5. output from running 19 gbk files (Design Set)  
`python3 moduleFile.py 19 1 boulderFile /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/gbks/ 19gbks_output1 20`  

It creates **3349** Ortho Groups, and for each, it generates 7 files:  
`OG0001903  primer3 config output`  
`OG0001903.fa  consambig output`  
`OG0001903.fasta multifasta file`  
`OG0001903Primers primer list (10 primers) for this OG`  
`OG0001903.ps primersearch result for the primers in this OG`  
`OG0001903.muscle  MUSCLE output`  
`OG0001903.trimAl  trimAl output`  

If I run `get_all_primers.py` to extract all the primer candidates it generated', the number count is 31774.  I'm sure there should be some filtering steps to remove a big chunk of these primers. But it's not clear to me what are those steps and how to run them.

**It seems that only around 2000 out of the 2461 primers can be found in this pool of 31774 primer candidates.**
  
***
**Mystery solved:**  
It turns out we have to use the older version of those dependencies in order to generate the primer pool which includes almost entirely those 2461 primers (or even those in the 7600 MasterPrimerFile), all but one **primer  GATATTCAAATTTGGCAGCCTGG AGATAGATCTGACTGGCGGTC   OG0001700primerGroup1**  

specifically, follow the steps below:  
1. Orthofinder 2.1.2  
download the binary file at this link: https://github.com/davidemms/OrthoFinder/tree/2.1.2#results-files  
And it's also required to install its own dependencies (module load them from scicomp):
    1. ml MCL-edge/12-068
    2. ml FastMe/2.1.5
    3. ml DLCpar/1.0
    4. ml ncbi-blast+/LATEST

2. TrimAl 1.2  
can't find where to download Trimal 1.2, so we use Trimal 1.4.1 instead of 1.2, you can create a conda env for this: `conda env create -n t3pio3 -f t3pio3.yaml`  

3. EMOBSS 6.4.0  
ml EMBOSS/6.4.0

4. Primer3 2.3.4  
ml primer3/2.3.4

5. MUSCLE 3.6  
we would use muscle/3.8.31 instead of MUSCLE 3.6.  (latest muscel version is 5.1)
I can't find where to download MUSCLE 3.6. (scicomp MUSCLE 3.6. doesn't work, missing other dependencies)  

6. once all the required packgaes are installed, run `python3 moduleFile_oldlibrary.py 19 1 boulderFile /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/gbks/ 19gbks_output1 20`  

**note:** changing the % of isolates inclusion will change the overal size of primer pool. But it seems that we can always find those 2460 primers in the final primer pool.  I tried this with 1 (all 19 isolates included), 0.93 (allowing 1 missing isolate) and 0.89 (allowing 2 missing isolates)   

***
### 6. filter part I  (run_primerscore.py)  

For the first pass filtering, I wrote a `run_primerscore.py` script, which utilizes Sean's `primerScore()` function. It selects one of the primer candidates (normally 10) in each oligo group, with the lowest score. If there is a tie, the first primer with the lowest score will be chosen.  The scripts seems to work well, except in the situation where there is a tie. `run_primerscore.py` always picks the first lowest score primer, but Sean either have additional filtering mechanism or he randomly picks one of the lowest score primers ?   

***
### 7. filter part II  (known genomes in Metagenomes and StoolBugs)   

In this step, we first use kraken2  

>  `for f in *.fasta;do kraken2 --db /scicomp/reference-pure/kraken_max $f  --threads 24 --use-names --output "${f%.fasta}.kraken";done`   

to classify all the contigs in those genomes and obtain a list of non-Salmonella contigs; While we run PrimerSearch `(run_primersearch.py)` to obtain the in silico amplicons of our primer candidates with those known genomes. We're then able to filter out `(parse_primersearch.py)` the primers that can amplify with non-Salmonella genomes, based on the list of non-Salmonella contigs and the PrimerSearch result. Below is the steps I tried:  
1. run `run_primerscore.py` as the first pass filtering and obtained a list of **3349** primers
2. run 2nd pass filtering as described above and obtained a list of **2822** primers.  
    - among these **2822** primers, **2235** of them match the **OG** part (Ex: **OG0001700**primerGroup9) of our **2461** primer panel.  
    - then among those **2235** primers, **1772** of them match completely, including the primer group number  

>  1.  I think if I reverse the filtering order, the result might be different. I didn't do that this time, because primersearch takes a long time (days) if we run it on a large number of primers (in our case, that will be **31774**).  
>  2.  Again, certainly I can try this and make the result closer to our current **2461** primer panel. But I'm also concerned the **2461** primer panel came out of running the pipeline with older version libraries. If we use all current version libraries, we will certainly get a different primer pool to start with.  
