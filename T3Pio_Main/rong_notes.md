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