# poolfileupload

This uploader assists in uploading Various RBTnSeq Files. 
* The Experiments files are TSVs with an undefined number of columns.
* The few required columns are: SetName, Index, Description, Date_pool_expt_started.
* For a full explanation of Experiments fcile refer to the end of the README.

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

---
* The biological goal of the program suite is to quantify the abundance of each 
barcoded strain in each sample.


* For more information, consider the following article:
"Rapid Quantification of Mutant Fitness in Diverse Bacteria by 
Sequencing Randomly Bar-Coded Transposons" - Kelly M. Wetmore, Morgan N. Price, 
[...], Adam P. Arkin, and Adam Deutschbauer
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4436071/







* Experiments Files:
* SetName -- experiments that were sequenced together should have the same set name, usually ending with set1, set2, etc.

* Index -- the multiplexing tag for sequencing, usually IT001, IT002, etc.

* Description -- Time0 means that this is a control sample, otherwise it describes the sample

* Date_pool_expt_started -- experiments and controls that were started on the same date are compared to each other.

Some other fields that are used if they are present:

* Condition_1 up to Condition_4 can be used to describe compounds that were added to the media. Units_1 and Concentration_1, etc., are used to indicate how much.

* Media can be used to describe the media. If it is used, then a warning is issued if it is not a known media (from feba/metadata/media)

* Group can be used to say which group of experiments this belongs to, i.e. "carbon source" or "nitrogen source". (Also, Group=Time0 means it is a control experiment.)

* The Media, Condition_*, and Group fields are included in some of the output tables, if they are present.


This is a [KBase](https://kbase.us) module generated by the [KBase Software Development Kit (SDK)](https://github.com/kbase/kb_sdk).

You will need to have the SDK installed to use this module. [Learn more about the SDK and how to use it](https://kbase.github.io/kb_sdk_docs/).

You can also learn more about the apps implemented in this module from its [catalog page](https://narrative.kbase.us/#catalog/modules/poolfileupload) or its [spec file]($module_name.spec).

# Setup and test

Add your KBase developer token to `test_local/test.cfg` and run the following:

```bash
$ make
$ kb-sdk test
```

After making any additional changes to this repo, run `kb-sdk test` again to verify that everything still works.

# Installation from another module

To use this code in another SDK module, call `kb-sdk install poolfileupload` in the other module's root directory.

# Help

You may find the answers to your questions in our [FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html) or [Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html).





