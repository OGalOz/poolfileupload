# poolfileupload

This uploader assists in uploading Various RBTnSeq Files. 
* The Experiments files are TSVs with an undefined number of columns.
* The few required columns are: SetName, Index, Description, Date_pool_expt_started.
* For a full explanation of Experiments fcile refer to the end of the README.
* Only upload multiple files at once if you're uploading PoolCount files.
* Create a single Output Name for each file you're uploading.
* If you're uploading any of 'ExperimentsTable', 'PoolCount' or 'PoolFile', then 
    you must indicate the GenesTable you are using. Otherwise, if you're uploading
    'InputGenesTable' - you do not use a reference GenesTable, since that is
    what you are creating.



---
* Pool Files are TSV files with 12 columns. 
* The headers are, in order from left to right: "barcode", "rcbarcode",
"nTot", "n", "scaffold", "strand", "pos", "n2", "scaffold2", "strand2", "pos2",
"nPastEnd".
* The pool file is part of the RBTnSeq (Random Barcode Transposon Sequencing)
analysis program. 
* The app "MapTnSeq" also outputs a pool file. 
* The pool file is used as a reference to find the important reads after 
running transposon sequencing. 
* A pool file is used as an input to the app "RBTnSeq Reads to Pool Counts".

---
* Pool Count Files are TSV files with a certain number of columns (TBD)
* The expected headers include: barcode, rcbarcode, scaffold, strand, pos 
* IT IS VERY IMPORTANT THAT THE OUTPUT NAME IS EXACTLY THE SET NAME FOR THE EXPERIMENT.
* Only import multiple TSV files if they are poolcount files!

---
* The biological goal of the program suite is to quantify the abundance of each 
barcoded strain in each sample.


* For more information, consider the following article:
"Rapid Quantification of Mutant Fitness in Diverse Bacteria by 
Sequencing Randomly Bar-Coded Transposons" - Kelly M. Wetmore, Morgan N. Price, 
[...], Adam P. Arkin, and Adam Deutschbauer
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4436071/






---
* Experiments Files:
Experiments files are tables that give information regarding the experiments. Each row
provides information about one experiment. So, for example, if you have 108 different
conditions tested (including the controls), then your Experiments table would have
108 rows not including the header row (with the column names). Column names
that are not optional:  "SetName", "Index", "Description", "Date_pool_expt_started", "Group"

* SetName -- experiments that were sequenced together should have the same set name, usually ending with set1, set2, etc. Note that the set and the index together should be unique for each experiment.

* Index -- the multiplexing tag for sequencing, usually IT001, IT002, IT025, etc. Note that the set and the index together should be unique for each experiment.

* Description -- "Time0" means that this is a control sample, otherwise it describes the sample - e.g. "N-Acetyl-D-Glucosamine nitrogen source".
                
* Date_pool_expt_started -- This should be a date value, e.g. "6/12/2021". Experiments and controls that were started on the same date are compared to each other.

* Group can be used to say which group of experiments this belongs to, i.e. "Time0", "carbon source" or "nitrogen source". (Also, Group=Time0 means it is a control experiment.)

Some other fields that are used if they are present:

* Drop -- This indicates whether or not to use an experiment. Values must be "True" or "False". "True" if you do not want to use the experiment.

* Condition_1 up to Condition_4 can be used to describe compounds that were added to the media. Units_1 and Concentration_1, etc., are used to indicate how much.

* Media can be used to describe the media. If it is used, then a warning is issued if it is not a known media (from feba/metadata/media)




There is an option to define whether an experiment is a "Time0" experiment (Control), 
    and which control group it belongs to by creating a column called "control_bool"
    where the values would be "True" if it is a control, and "False" otherwise. 
    And then you would also need to create a column called "control_group"
    where the experiments which
---
    


This is a [KBase](https://kbase.us) module generated by the [KBase Software Development Kit (SDK)](https://github.com/kbase/kb_sdk).


You can also learn more about the apps implemented in this module from its [catalog page](https://narrative.kbase.us/#catalog/modules/poolfileupload) or its [spec file]($module_name.spec).



