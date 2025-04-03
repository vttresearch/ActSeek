![image](https://github.com/user-attachments/assets/92532fa4-0cc2-4f10-9c2f-f84e76a137a9)# ActSeek

ActSeek is a tool for enzyme mining in the Alphafold database. ActSeek uses the coordinates of at least 3 user defined amino acids from a seed structure to find other proteins with the same (or similar) amino acids in the same position. 

## Installation

### Prerequisites

- Python 3.6 or higher
- `pip` (Python package installer)

### Required Libraries

Install the required libraries using `pip`:

```sh
pip install numpy biopython tqdm wget requests
pip install cython
pip install scipy
pip install pyKVFinder
pip install .

```
The packages required are:

numpy
biopython
tqdm
wget
requests
scipy
pyKVFinder
cython

## Usage

To test ActSeek after its installation:
```sh
cd tests
actseek

```

It will create a folder "results" with the results of the search. 

Doing you own search:
Edit the json configuration file "config.json". 

### comman line arguments

The program accepts several command line arguments to override the configuration file:

-a, --active-site: Position of the amino acids in the seed structure involved in the search.

-sa, --selected-active: Selection of the 3 main amino acids involved in the search. Default: 0,1,2

-g, --aa-grouping: Amino acid grouping. Dictionary of amino acids and their correspoing group. the vaule written as the group is not relevant but all the amino acids belonging to the same group should have the same value. 

-r, --random-seed: Random seed.

-t1, --threshold: Threshold of the average distance of the mapped amino acids.

-t1c, --threshold-combinations: Threshold for choosing the possible hits. Usually bigger than the t1 threshold.

-t2, --aa-surrounding: Number of amino acids around the selected amino acids that will be taken into account while deciding which is the active site.

-t2t, --aa-surrounding-threshold: Maximum average distance of the surrounding amino acids to their respective matching amino acids to be selected.

-t3, --threshold-others: Maximum distance to the matching amino acid to be considered in the final solution.

-i, --iterations: Maximum number of iterations.

-f, --first-in-file: First protein to perform tests with in the file (Usefull when it is run in batch mode).

-m, --max-protein: Max number of proteins to perform tests with (Usefull when it is run in batch mode).

-s, --protein-file: Path to the file containing the name of the proteins to be tested.

-af, --alphafold-proteins-path: Path to Alphafold proteins folder.

-p, --seed-protein-file: Path to seed protein file.

-d, --delete-protein-files: Delete protein files.

-pr, --path-results: Path of the results.

-ts, --testing: Testing one protein. This argument takes the Uniprot Id of the protein and returns the structure of the protein aligned with the seed structure. 

-kv, --KVFinder:Uses KVFinder to compare the cavity where the active sides are with the cavity in the seed structure.

-c, --custom: Using a custom structure database not connected with Alphafold for the search.

#### Comman line example
```sh
actseek -a "292_A,448_A,478_A" -s test.txt -p AF-Q9ZHI2-F1-model_v4.pdb -pr "results" -t1 3 -kv
```
I you are using more than 3 amino acids, remember to add the "selected_active" parameter:
```sh
actseek -a "145_A,180_A,292_A,448_A,449_A,478_A" -sa "2,3,5" -s test.txt -p AF-Q9ZHI2-F1-model_v4.pdb -pr "results" -t1 3 -kv
```
#### Testing one protein
If you are interested on getting the output mapping file for one specific protein, the option -ts can be used followed by the id of the protein you want to test. For example:
```sh
actseek -a "292_A,448_A,478_A" -s test.txt -p AF-Q9ZHI2-F1-model_v4.pdb -pr "results" -t1 3 -ts P45370
```
This will output the mapping file and the pdb file of the protein P45370 superposed with the seed protein (in this case Q9ZHI2). The alignment is based on the ActSeek algorithm using only the selected amino acids as mapping seed. 

![image](https://github.com/user-attachments/assets/1c94d140-a78f-45ac-8630-69bd26ec38d0)
To visualize the protein P45379 superpostion with the seed protein, I open the file P45279.pdb created by ActSeek in the results folder and the seed structure using Pymol.

### Needed files
- Seed pdb file (in tests it is AF-Q9ZHI2-F1-model_v4.pdb)
- List of Uniprot ids (in tests it is test.txt)
- config.json file for the configuration that can be introduced also in the command line.

### Config.json

{

    "active_site": "292_A,448_A,478_A",
    
    "selected_active": "0,1,2",
    
    "aa_grouping": {
        "GLY": "GLY",
        "ALA": "ALA",
        "PRO": "PRO",
        "ARG": "ARG",
        "HIS": "HIS",
        "LYS": "LYS",
        "ASP": "ASP",
        "GLU": "GLU",
        "SER": "SER",
        "THR": "THR",
        "ASN": "ASN",
        "GLN": "GLN",
        "CYS": "CYS",
        "VAL": "VAL",
        "ILE": "ILE",
        "LEU": "LEU",
        "MET": "MET",
        "PHE": "PHE",
        "TYR": "TYR",
        "TRP": "TRP",
        "HSD": "HIS"
    },
    
    "random_seed": 0,
    
    "threshold": 1.0,
    
    "threshold_combinations": 3.0,
    
    "aa_surrounding": 4,
    
    "aa_surrounding_threshold": 3.0,
    
    "threshold_others": 4.0,
    
    "iterations": 2000,
    
    "first_in_file": 0,
    
    "max_protein": 10000,
    
    "protein_file": "test.txt",
    
    "alphafold_proteins_path": "structures",
    
    "seed_protein_file": "AF-Q9ZHI2-F1-model_v4.pdb",
    
    "delete_protein_files": true,
    
    "path_results": "results",
    
    "testing": "",

    "KVFinder": false,

    "custom": false
    
}

To perform a search where two or more amino acids can be in interchangeble, change the values of the aa_grouping variable to be the same.

### Using a customized database for the search

All the structures from the customized set can be downloaded into a folder that will be given to ActSeek using the parameter "alphafold_proteins_path". 
The "protein_file" should contain the name of the proteins that you want to use in the search without the extension (one name in each line). 
run actseek with the parameter -c:
```
actseek -c -af "your_protein_path" -s "your_protein_list"
```
