# DiMER
DiMER is a modular protein annotation pipeline that allows the user to query their protein sequences against any number of BLAST and/or HMM databases. Notably, DiMER also allows the user to input any number of case-insensitive text strings they’d like to be preferentially avoided in the annotation of their submitted proteins (ex: "hypothetical protein").

## Dependencies

- [Python](https://www.python.org/) (Tested with version 3.13.1)
- [HMMER3](http://hmmer.org/) (Tested with version 3.1b2)
- [Diamond](https://github.com/bbuchfink/diamond) (Tested with version 2.1.10)
- Python Library: [Pandas](https://pandas.pydata.org/)

- Note, the workflow of DiMER requires editing some files in this program. To do this, we recomend using a plain text editor. You can certainly use vim or nano from your terminal, but if you're not as comfortable using the command line, you may enjoy using an IDE like [VScode] instead.

## Getting Started

### Installing Dependencies
The standalone dependencies of DiMER can be installed as system-wide executables and DiMER *should* work. However, because HMMER3 and Diamond have their own dependencies, we recomend installing DiMER dependencies in a virtual environment using [conda](https://docs.anaconda.com/miniconda/) as it will handle acquiring and managing all software dependencies for you.

__Note__: you can use newer versions of Diamond and HMMER assuming you are using databases formated to those newer versions (see formatting instructions below in __Setting up your Databases__). If you are using *our* MGE database collection you will need HMMER 3.1b2 and Diamond 2.2.10.

1. Using Conda, you will create an environment for DiMER and activate it before each use of the program.
- Make sure you have the channels bioconda and conda-forge added to your conda config file
```
conda config --add channels bioconda
conda config --add channels conda-forge
```
- Create an environment for DiMER with Python 3.13.1 included.
```
conda create -n dimer python=3.13.1
```
- Activate the DiMER environment. You will run this command each time before using DiMER.
```
conda activate dimer
```

2. Once you've activate the conda environment, you can begin installing the dependencies. You will only need to do this once.


- install HMMER3
```
conda install bioconda::hmmer=3.1b2
```
- install Diamond
```
conda install diamond=2.1.10
```
- install the Python library Pandas
```
conda install pandas
```

__Note:__ If you're new to conda, be aware of error messages during install. These errors are specific to your machine and may need to be addressed before these packages will install and work correctly.

3. Clone the DiMER project repository to a directory of your choosing.
- Download DiMER (requires git (conda install anaconda::git) or you can download from your browser here: https://github.com/sandialabs/DiMER/tree/main)
```
git clone https://github.com/sandialabs/DiMER.git
```
### Setting up your Databases
DiMER was developed to annotate the proteins from mobile genetic elements in prokaryotes. However, it can be used for any taxa assuming the appropriate databases are created/downloaded and formatted. 

__Installing MGE Databases__
1. If you would like to use the Diamond (2.1.10) and HMMER (3.1b2) formatted databases tailored to MGE annotation that we used in our publication you can download them from your browser at https://zenodo.org/records/15724972 and then move them to your dimer directory. Or, use the included python script to fetch and untar the databases and associated files:

- From within the dimer directory run (you will need ~25GB of storage available):
```
python get_databases.py
```

__Note 1__: If downloading from a browser remember to extract the tar archives ("tar -xzvf <file/folder_name>.tar.gz") after download to your dimer folder.

__Note 2:__ If you decide you *don't* like using one or more of these databases you can add a '#' to the beginning of the line containing the path to the database in either 'hmmer3_db_paths.txt' or 'diamond_db_paths.txt' to exclude the database from being queried in your run.

__Using your Own Databases__
You can use your own collection of protein databases provided you format them for use with Diamond and HMMER3.
1. To generate a BLASTP database to use with Diamond, you will need all proteins for the database in a single file in FASTA format. Then, to generate the database you'll run:

```
diamond makedb --in <fasta file> -d <name of your new database>
```
2. To generate a HMM for use with HMMER you'll use a multiple sequence alignment of a single protein to create a single HMMER formatted HMM profile:
```
hmmbuild <msa_file.out> <msa file>
```
- However, you probably want more than one hmm profile in your hmm database. Assuming all of your hmm profiles have the same subcript, you can concatenate them all into a single hmm database using cat:
```
cat *.out > <hmm database_name>.hmm
```
3. Once your databases are made you'll need to add their paths to diamond_db_paths.txt or hmmer3_db_paths.txt. Make sure to follow the same format. Alternatively, you can create your own files with the same format and give them to dimer with the flag ```-blast_paths <path to file>``` and/or ```-hmm_paths <path to file>```
-An example of the diamond_db_paths.txt format is below:
```VFDB=<full_path_to_database_file>/VFDB_setB_pro.dmnd```
-An example of the hmmer3_db_paths.txt format is below:
```Pfam37A=<full_path_to_database_file>/Pfam-A.hmm```
Here, your desired database name preceeds the equal sign. This database name will be carried over to DiMER outputs and it is important this name exactly matches the first column entries in the Database_key.tsv file if you need to use it (see next section).

4. Some databases use a key:value table to provide annotations. An example is the [PHROG](https://phrogs.lmge.uca.fr/READMORE.php) database that annotates their database with values such as PHROG1 and PHROG123 which correspond to the true annotation values (integrase and tail-terminator connector) in another table. If you would like these keys automatically replaced with their respective values in the final output of dimer you'll need to populate the file Databases_key.tsv which you can do in excel, if you're not comfortable with coding, as long as you save it as a tsv. There are four columns in this table. The first column correspnds to the name you've given the database in the path file. The second is the key given by the annotation table from your database and the third is the value that will be carried over in the final annotation of your protein (if it's a match, of course). The last column is optional to *populate* (though it must still exist!) and will be carried over in the DiMER output table containing all annotation matches (combined_output.tsv) but not your final re-annotated protein file. This is a good place to put broad functional classifications if they exist for the database you're using and they might be of interest to your research.
Here's an example of the populated table:
```
PHROG	phrog_1	integrase	integration and excision
PHROG	phrog_123	tail terminator	connector
CAZy	AA8	iron reductase domain	Auxiliary activities
CAZy	AA9	copper-dependent lytic polysaccharide monooxygenases 	Auxiliary activities
CAZy	CBM0	unclassified	Carb binding modules
CAZy	CBM10	cellulose binding	Carb binding modules
```

Once you've populated this file you'll need to run the included script to turn this table into a "pickle" for DiMER to use:
```
python repickle.py
```
 The table must retain the name "Databases_key.tsv" for this script to work.
### Editing Placeholder Annotations
The "special" thing that sets DiMER apart from other annotation tools is that it attempts to selectively exclude placeholder annotations (*e.g.:* "hypothetical protein", “uncharacterized protein”, “unannotated protein”, “conserved ORF”, “predicted protein”, “putative ORF”, “unnamed protein product”, “unknown function”, “ORF”-followed by a series of integers, “putative coding sequence”, or “hypothetical ORF” and some taxon-specific placeholders —such as “putative phage protein” or “hypothetical bacterial protein”) assuming there is a match within your databases with a sufficient e-value to replace them. To include or exclude terms to be selectively excluded from annotation with DiMER you'll need to edit the file "placeholders.txt". At download, this file contains some common placeholder terms and you can continue adding to this file or make your own file and include it in your DiMER run with the flag ```-placeholer <path to placeholder file>```. Notably, by including the character '#' DiMER will be able to replace any set of integers for terms that might include numbers. For example, adding ORF# to the list will be able to replace instances of ORF11, ORF4056, ORF39, etc. The entries in placeholders.txt are not case sensitive, but they are sensitive to word order, punctuation, and spacing.

___Warning___: Anything you include in placeholder.txt will be excluded so it's important that you're considerate and careful when making this list. For example, if I chose to use only the work "hypothetical" I would exclude hypothetical protein in my final output but also potentially useful annotations like "hypothetical RNA polymerase subunit". Thus it is also important you leave no blank lines or lines containing only a space character.

## Running DiMER
### Checklist:
__Before running DiMER make sure you have:__ 
1) Activated your conda environment: ```conda activate dimer```.
2) Either downloaded our pre-formated [MGE_Databases](https://zenodo.org/records/15724972) or installed and formatted your own databases with Diamond and HMMER.
3) Listed the paths to the databases you'd like to query in either hmmer3_db_path.txt and diamond_db_paths.txt (if you're using our MGE_databases this is already done and included) or another set of files in the same format.
4) Generated Databases_key.pickle if applicable to your dataset (if you're using our MGE_databases this is already done and included).
5) Included any placeholder annotations such as "hypothetical protein" that you'd like to be selectively avoided in the annotations of your proteins in the file placeholders.txt. __Remember:__ leave no empty lines in this file. Empty annotations and instances of "NA" are already handled by the program.
5) Obtained a set of translated proteins in fasta format that you'd like to annoate, or re-annotate, using DiMER!

### Quickstart
Requirements are an input file containing translated coding sequences to be annotated in fasta format and an output directory for your results to be placed in.
```
python dimer.py -in <path to translated proteins in fasta format> -out <output directory>
```
DiMER will create the output directory you specify if it does not exist. __Note:__ it is recomended to use a unique output directory as the -restart flag will delete all contents of this directory if used.

Running DiMER this was will only create the tsv table of all the best hit annotations for each protein across each database. This file lists all best annotations but not necessarily *good* annotations. Pay attention to the e-values and ensure they're at least low enough to be considered credible!

### Recommended Run
 If you'd like to generate a fasta file with your proteins where you selectively exclude placeholder annotations and only use annotations below a specified e-value (Default<=10^-10 or change this value with the ```-eval``` flag) you will need to use the ```-replace``` flag with the option ```force``` or ```gentle```.

 ___1. If your input protein file is already annotated and you just want to replace placeholder annotations with better options___ (if they exist and are below the evalue threshhold): use ```-replace gentle``` to only replace placeholder annotations present in your protein fasta file (you can add or subtract from the list of placeholder annotations by editing the file placeholder.txt).
```
python dimer.py -in <path to translated proteins in fasta format> -out <output directory> -replace gentle
```
A new fasta file will be generated in your DiMER output folder called ```<your fasta file name>.updated.faa``` unless you specify a new name with the flag ```-unique_name```.

___2. If your input file is *not* annotated or if you would like to replace *all* existing annotations with DiMER___ (assuming they exist in your chosen databases and are below the evalue threshhold): use ```-replace force``` to generate new annotations or replace *all* annotations present in your original file. This mode will also use placeholders.txt to selectively exclude any placeholder annotations listed therein whenever possible.

```
python dimer.py -in <path to translated proteins in fasta format> -out <output directory> -replace force
```
Again, a new fasta file will be generated in your DiMER output folder called ```<your fasta file name>.updated.faa``` unless you specify a new name with the flag ```-unique_name```.

### DiMER Flags and Options:


- ```-clean```


    All created folders within the specified output directory will be deleted. This will remove extraneous BLAST and HMM results to free up sapce on your computer. I would only do this if you're not planning to have to add any placeholder annotations and rerun dimer since the database searching won't be repeated if that databases output folder is present in the DiMER output folder.

- ```-restart```

    All created folders within the specified output directory will be deleted as a prep to freshly restart the run. This is suggested if the run originally fails due to a terminal session close or compute/RAM restrictions.

- ```-unique_name <new output file name>```

    Specify a file name instead of the default "combined_output.tsv"

- ```-best_hits```

    Outputs a file called best_hits.txt (unless you specify a unique name then it will be "<unique_name>".best_hits.txt) which contains a list of the annotations with the lowest e-value score with preference given to low scoring hits that are not listed in the placeholders.txt file.

- ```-replace <force or gentle>```

    - ```-replace gentle```: Replace hypothetical annotations in the original translated protein fasta file file with new annotations. "Gentle" will only replace if the annoation string contains a phrase from placeholder.txt. Do *not* use this mode on files that have not been previously annotated. 
    - ```-replace force```: Reannotates all sequences in your translated fasta file with the best scoring hit across all databases queried while selectively avoiding entries that match those in placeholders.txt whenever possible. 
    - Both options will create a new file with your original input file prefix followed by ".updated.faa" in the specified output dir.

- ```-eval <specify maximum evalue threshold>```

    - Specify the maximum evalue threshold for inclusion in the list of best hits for final annotation of prot file in "-replace". Using this flag without using the "-replace" or "-best_hits" flag will have no effect. Default: 1e-10

- ```-cpu <number of cpus>```

    - Specify the number of additional cpus to use for HMMER - the default is zero. Must be an integer.

- ```-placeholders <path to file containing list of placeholder annotations>```

    - Specify a file of placeholder annotations to selectively replace in your annotation run. This might be useful if you use DiMER for different projects or on different taxa with different common placeholder annotations. Only one entry should be made per each line in the file. No line should be empty. Default runs use "placeholders.txt". We recommend using this file as a template.


## Citations

If using this software please cite:

__The DiMER Publication__

PAPER CITATION WHEN WE HAVE IT

__Diamond__
- Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

__Hmmer3__
- HMMER 3.4 (Aug 2023); http://hmmer.org/
Copyright (C) 2023 Howard Hughes Medical Institute.
Freely distributed under the BSD open source license.

If using our formated MGE_Databases please cite the publications they come from:

- Terzian, P., et al., PHROG: families of prokaryotic virus proteins clustered using remote homology. NAR Genom Bioinform, 2021. 3(3): p. lqab067.
- Brown, C.L., et al., mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. Appl Environ Microbiol, 2022. 88(18): p. e0099122.
- Wang, M., et al., ICEberg 3.0: functional categorization and analysis of the integrative and conjugative elements in bacteria. Nucleic Acids Res, 2024. 52(D1): p. D732-D737.
- Payne, L.J., et al., Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. Nucleic Acids Res, 2021. 49(19): p. 10868-10878.
- Mistry, J., et al., Pfam: The protein families database in 2021. Nucleic Acids Res, 2021. 49(D1): p. D412-D419.
- Cury, J., et al., Identifying Conjugative Plasmids and Integrative Conjugative Elements with CONJscan. Methods Mol Biol, 2020. 2075: p. 265-283.
- Feldgarden, M., et al., AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep, 2021. 11(1): p. 12728.
- Trgovec-Greif, L., et al., VOGDB-Database of Virus Orthologous Groups. Viruses, 2024. 16(8).
- Haft, D.H., J.D. Selengut, and O. White, The TIGRFAMs database of protein families. Nucleic Acids Res, 2003. 31(1): p. 371-3.


## Copywrite
Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
