# ActSeek

ActSeek is a tool for enzyme mining in the Alphafold database. ActSeek uses the coordinates of at least 3 user defined amino acids from a seed structure to find other proteins with the same (or similar) amino acids in the same position. 

## Installation

### Prerequisites

- Python 3.6 or higher
- `pip` (Python package installer)

### Required Libraries

Install the required libraries using `pip`:

```sh
pip install .

```
The packages required are:

numpy

biopython

tqdm

wget

requests

scipy

cython

### Usage

To test ActSeek after its installation:
```sh
cd tests
actseek

```

It will create a folder "results" with the results of the search. 

Doing you own search:
Edit the json configuration file "config.json". 

#### comman line arguments

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


