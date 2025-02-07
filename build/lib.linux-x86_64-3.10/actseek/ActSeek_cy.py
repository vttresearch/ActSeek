import numpy as np
import random
from Bio.PDB import *
from Bio.PDB import PDBParser
import warnings
import os
import argparse
import concurrent.futures
from multiprocessing import Process
import actseek.ActSeekLib_cy as ActSeekLib
import traceback
from tqdm import tqdm
import wget
import requests
import json
import tempfile
import unittest

warnings.filterwarnings("ignore")

config = None


# Reads the protein structure with the Biopython library and extracts the coordinates of the amino acids.
def read_pdbs(aa_active, seed_protein, case_protein):
    """
    Reads PDB files for the seed and case proteins, extracts coordinates and other relevant information.

    Parameters:
    aa_active (list): List of active amino acids in the seed protein.
    seed_protein (str): Path to the seed protein PDB file.
    case_protein (str): Path to the case protein PDB file.

    Returns:
    tuple: A tuple containing the following elements:
        - protein_coords (np.ndarray): Coordinates of the alpha carbon atoms (CA) in the case protein.
        - protein_coords_cb (np.ndarray): Coordinates of the beta carbon atoms (CB) in the case protein.
        - seed_coords (np.ndarray): Coordinates of the alpha carbon atoms (CA) in the seed protein.
        - cavity_coords (np.ndarray): Coordinates of the alpha carbon atoms (CA) in the active site of the seed protein.
        - cavity_coords_cb (np.ndarray): Coordinates of the beta carbon atoms (CB) in the active site of the seed protein.
        - aa (dict): Dictionary mapping residue indices to residue names in the case protein.
        - aaCav (dict): Dictionary mapping residue indices to residue names in the active site of the seed protein.
        - active (np.ndarray): Array of active residue indices in the seed protein.
        - real_index (dict): Dictionary mapping sequential indices to residue indices in the case protein.
        - real_index_seed (dict): Dictionary mapping residue indices to sequential indices in the seed protein.    
    """
    parser = PDBParser()
    seed_structure = parser.get_structure("complex", seed_protein)
    case_structure = parser.get_structure("complex2", case_protein)
    aa = {}
    protein_coords = []
    protein_coords_cb = []
    i = 0
    real_index = {}
    for model in case_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if "CA" in atom.fullname:
                        protein_coords.append(atom.get_coord())
                    if 'CB' in atom.fullname:
                        protein_coords_cb.append(atom.get_coord())
                if residue.get_resname() == "GLY":
                    protein_coords_cb.append([-10000000, -10000000, -10000000])
                aa[int(residue.get_id()[1])] = str(residue.get_resname())
                real_index[i] = int(residue.get_id()[1])
                i = i + 1
            break
    cavity_coords = []
    cavity_coords_cb = []
    seed_coords = []
    aaCav = {}
    active = []
    real_index_seed = {}
    i = 0
    for model in seed_structure:
        for chain in model:
            for residue in chain:
                real_index_seed[int(residue.get_id()[1])] = i
                i = i + 1
                for atom in residue:
                    if "CA" in atom.fullname:
                        seed_coords.append(atom.get_coord())

                if residue.get_id()[1] in aa_active:
                    active.append(int(residue.get_id()[1]))
                    for atom in residue:
                        if "CA" in atom.fullname:
                            cavity_coords.append(atom.get_coord())
                        if 'CB' in atom.fullname:
                            cavity_coords_cb.append(atom.get_coord())
                    if residue.get_resname() == "GLY":
                        cavity_coords_cb.append([-10000000, -10000000, -10000000])
                    aaCav[int(residue.get_id()[1])] = str(residue.get_resname())

    
    return np.array(protein_coords), np.array(protein_coords_cb), np.array(seed_coords), np.array(cavity_coords), np.array(cavity_coords_cb), aa, aaCav, np.array(active), real_index, real_index_seed



def printProtein(translation_vector, rotation, path, name):
    """
    Transforms and saves a protein structure by applying a given rotation and translation.

    Args:
        translation_vector (numpy.ndarray): A 3-element array representing the translation vector.
        rotation (numpy.ndarray): A 3x3 matrix representing the rotation to be applied.
        path (str): The file path to the input PDB file.
        name (str): The name to be used for the output PDB file.

    Returns:
        None
    """
    parser = PDBParser()
    structure2 =  parser.get_structure("complex2", path)

    for model in structure2:
        for chain in model:
            for residue in chain:
                for atom in residue:
                     atom.transform(rotation, translation_vector)

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure2)
    io.save(config.path_results+"/"+name+".pdb")

def ActSeekMain(aa_active, aa_des, seed_protein, case_protein_filename, iterations, case_protein_name):
    """
    Main function for the ActSeek algorithm.

    Parameters:
    aa_active (list): List of active amino acids.
    aa_des (dict): Dictionary of desired amino acids.
    seed_protein (str): Filename of the seed protein structure.
    case_protein_filename (str): Filename of the case protein structure.
    iterations (int): Number of iterations for the RANSAC algorithm.
    case_protein_name (str): Name of the case protein.

    Returns:
    tuple: Contains the following elements:
        - finalDist (float): Final distance after RANSAC.
        - distances (list): List of distances for the solution.
        - distances_arround (float): Distance around the solution.
        - t_transformed (numpy.ndarray): Transformed coordinates.
        - solution (list): List of valid amino acid correspondences.
        - translation_vector (numpy.ndarray): Translation vector for the transformation.
        - rotation (numpy.ndarray): Rotation matrix for the transformation.
    """
    # Reads the structures of the seed and the case structures with the Biopython library and extract the needed information

    try:
        protein_coords, protein_coords_cb, seed_coords,cavity_coords, cavity_coords_cb,  aa, aaCav, active, real_index, real_index_seed = read_pdbs(
            aa_active,
            seed_protein,
            case_protein_filename)
    except:
        traceback.print_exc()
        return None, None, None, None, None, None, None, None, None, None

    index_used = np.array([int(x) for x in config.selected_active.split(",")])
    active_used = np.array([active[x] for x in index_used])
    cavity_coords_used = np.array([cavity_coords[x] for x in index_used])
    cavity_coords_cb_used = np.array([cavity_coords_cb[x] for x in index_used])

    # Creates an object of the class Active_side where all the information is added
    active_site = ActSeekLib.Active_site(active_used, aaCav, cavity_coords_used, aa_des)
    active_site.get_AA_groups(False)

    # Produces a vector with all the possible amino acid correspondences based on the aa_des dictionary
    pc = active_site.get_possible_correspondences(protein_coords, aa, real_index)
    threshold = config.threshold_others

    # Calculate the valid combinations by removing all the combinations of 3 amino acid mapping that are not in the appropiated distances with a threshold of 1 (that coud be an argument for the main function to be defined by the user)
    valid_combinations = active_site.get_all_possible_combinations(pc, protein_coords, aa, real_index, threshold,
                                                                   aa_des)


    if len(valid_combinations) > 0:
        distances, distances_arround, t_transformed, solution, translation_vector, rotation = ActSeekLib.ransac_protein(
            cavity_coords, protein_coords, protein_coords_cb, seed_coords, valid_combinations,
            iterations, aa_des, aa, real_index, real_index_seed, aaCav, active, cavity_coords_used,cavity_coords_cb_used, config.aa_surrounding, config.threshold_combinations, 0)


        if np.sum(distances) / len(solution) < config.threshold and len(solution) >= 3 and distances_arround < config.aa_surrounding_threshold:
            try:
                rmsd, minrmsd, percentage= ActSeekLib.getGlobalDistance(t_transformed, seed_coords)
            except:
                traceback.print_exc()

            if len(config.testing) > 4:
                sol_write = open(f"{config.path_results}/{case_protein_name}.txt", "w")
                printProtein(translation_vector, rotation, case_protein_filename, case_protein_name)
            else:
                sol_write = open(f"{config.random_dir}/{case_protein_name}.txt", "w")

            strmapping = ""
            for aamap in solution:
                aa_protein = aa[real_index[aamap[0]]]
                aa_cavity = aaCav[active[aamap[1]]]
                strmapping = strmapping + str(real_index[aamap[0]]) + aa_protein + "-" + str(
                    active[aamap[1]]) + aa_cavity + ";"


            sol_write.write(
                case_protein_name + "," + strmapping + "," + str(np.sum(distances) / len(active)) + "," + str(distances_arround) + "," + ";".join(map(str, distances)) + ","+str(rmsd)+","+str(minrmsd)+","+str(percentage)+"\n")
            sol_write.close()
    else:
        pass



def processProtein(case_protein_name):
    seed_protein = config.seed_protein_file
    case_protein_name = case_protein_name

    try:
        # Download the protein structure from alphafold database
        case_protein_filename = f"AF-{case_protein_name}-F1-model_v4.pdb"
        case_protein_filepath = f"{config.alphafold_proteins_path}/{case_protein_filename}"

        if os.path.isfile(case_protein_filepath) == False:
            try:
                url= f"https://alphafold.ebi.ac.uk/files/{case_protein_filename}"
                response = requests.head(url)
                response.raise_for_status()
                wget.download(url, out=case_protein_filepath, bar=None)
            except requests.exceptions.RequestException:
                pass

        # This part is give by the user, it should maybe be in the main function to be given as argument
        aa_active = [int(x) for x in config.active_site.split(",")]
        

        # Main function of the algorithm
        ActSeekMain(aa_active, config.aa_grouping, seed_protein, case_protein_filepath, config.iterations, case_protein_name)

        # Removes the protein structure once the algorithm has finished
        if config.delete_protein_files == True:
            try:
                os.remove(case_protein_filepath)
            except:
                traceback.print_exc()

                pass
    except:
        traceback.print_exc()


def read_config(file_path):
    with open(file_path, 'r') as file:
        config = json.load(file)
    return config




def parse_args():
    config = read_config('config.json')    

    parser = argparse.ArgumentParser(description="ActSeek Program")
    parser.add_argument("-a", "--active-site", type=str, default=config['active_site'],
                        help="Position of the amino acids involved in the search")
    parser.add_argument("-sa", "--selected-active", type=str, default=config['selected_active'],
                        help="Position of the three amino acids involved in the search.")
    parser.add_argument("-g", "--aa-grouping", type=json.loads, default=json.dumps(config['aa_grouping']),
                        help="Amino acid grouping")
    parser.add_argument("-r", "--random-seed", type=int, default=config['random_seed'], help="Random seed")
    parser.add_argument("-t1", "--threshold", type=float, default=config['threshold'],
                        help="Threshold of the average distance of the mapped amino acids.")
    parser.add_argument("-t1c", "--threshold-combinations", type=float, default=config['threshold_combinations'],
                        help="Threshold for choosing the possible hits.")
    parser.add_argument("-t2", "--aa-surrounding", type=int, default=config['aa_surrounding'],
                        help="Number of amino acids around the selected amino acids that will be taken into account while deciding which is the active site.")
    parser.add_argument("-t2t", "--aa-surrounding-threshold", type=float, default=config['aa_surrounding_threshold'],
                        help="Maximum average distance of the surrounding amino acids to their respective matching amino acids to be selected.")
    parser.add_argument("-t3", "--threshold-others", type=float, default=config['threshold_others'],
                        help="Maximum distance to the matching amino acid to be considered in the final solution.")
    parser.add_argument("-i", "--iterations", type=int, default=config['iterations'], help="Maximum number of iterations")
    parser.add_argument("-f", "--first-in-file", type=int, default=config['first_in_file'], help="First protein to perform tests with")
    parser.add_argument("-m", "--max-protein", type=int, default=config['max_protein'],
                        help="Max number of proteins to perform tests with")
    parser.add_argument("-s", "--protein-file", type=str, default=config['protein_file'],
                        help="Path to the file containing the name of the proteins to be tested")
    parser.add_argument("-af", "--alphafold-proteins-path", type=str, default=config['alphafold_proteins_path'],
                        help="Path to Alphafold proteins folder")
    parser.add_argument("-p", "--seed-protein-file", type=str, default=config['seed_protein_file'],
                        help="Path to seed protein file")
    parser.add_argument("-d", "--delete-protein-files", dest="delete_protein_files", action="store_true", default=config['delete_protein_files'],
                        help="Delete protein files")
    parser.add_argument("-pr", "--path-results", type=str, default=config['path_results'], help="Path of the results")    
    parser.add_argument("-ts", "--testing", type=str, default=config['testing'],
                    help="Testing one protein. This argument takes the Uniprot Id of the protein.")
    
    args = parser.parse_args()
   
    return vars(args)


def main():
    # Get experiment config from command line arguments
    global config
    config = parse_args()
    if os.path.isdir(config['path_results']) == False:
        os.mkdir(config['path_results'])

    if os.path.isdir(config['alphafold_proteins_path']) == False:
        os.mkdir(config['alphafold_proteins_path'])

    with tempfile.TemporaryDirectory(dir=config['path_results']) as tmpdir:
        config['random_dir'] = tmpdir
        

        config = argparse.Namespace(**config)

        # Set random seed to a fix value, to ensure repeatability
        random.seed(config.random_seed)

       
        if len(config.testing) > 4:
            case_protein = config.testing
            p = Process(target=processProtein, args=(case_protein,))
            p.start()
            p.join()
        else:
            # The cluster will send jobs with a number from 0 to 250 (2500000/100000) and each job will process 10000 protein
            init = int(config.first_in_file) * 10000

            file = open(config.protein_file, "r")

            start = False
            cases =[]
            for i, line in enumerate(file):
                if i == init:
                    start = True
                if start:
                    case_protein = line.split('\t')[0]
                    case_protein = case_protein.split('\n')[0]
                    cases.append(case_protein)
                if i == init + config.max_protein+2:
                    break
            file.close()

            
            with tqdm(total=len(cases)) as pbar:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    #future_to_item = {executor.submit(processProtein, case): case for case in cases}
                    futures = [executor.submit(processProtein, case) for case in cases]
                    for future  in concurrent.futures.as_completed(futures):
                        pbar.update(1)

            files = os.listdir(config.random_dir)
            results = open(config.path_results+"/results.csv","w")
            results.write("Uniprot ID,Mapping,Average distance,Average distance AA arround, All distances,Structural similarity, Structural RMSD, Percentage structural mapping\n")
            for file in files:
                f = open(config.random_dir+"/"+file, "r")
                for line in f:
                    results.write(line)
                f.close()
                os.remove(config.random_dir+"/"+file)
            results.close()


if __name__ == '__main__':
    main()
