import numpy as np
import random
from Bio.PDB import PDBParser, PDBIO
import warnings
import os
import argparse
import concurrent.futures
import actseek.ActSeekLib_cy as ActSeekLib
import traceback
from tqdm import tqdm
import wget
import requests
import json
import tempfile
from functools import partial
import pyKVFinder

warnings.filterwarnings("ignore")

config = None


def read_pdbs_case(case_protein_path):
    """
    Reads a PDB structure file and extracts alpha carbon (CA) and beta carbon (CB) coordinates,
    while maintaining residue names and an index mapping.

    Args:
        case_protein_path (str): The file path to the PDB structure.

    Returns:
        tuple:
            - np.ndarray: Coordinates of CA atoms for all residues.
            - np.ndarray: Coordinates of CB atoms for all residues (with a placeholder for Gly).
            - dict: A mapping from residue index (as in the PDB file) to residue name.
            - dict: A global index mapping from the function’s residue order to the PDB’s residue index.

    Note:
        For Glycine (GLY), a placeholder coordinate [-10000000, -10000000, -10000000] 
        is used for the missing CB atom.
    """
    
    pdb_parser = PDBParser()
    pdb_structure = pdb_parser.get_structure("complex2", case_protein_path)

    residue_name_map = {}
    ca_coords = []
    cb_coords = []
    residue_counter = 0
    global_res_index = {}

    for model in pdb_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if "CA" in atom.fullname:
                        ca_coords.append(atom.get_coord())
                    if "CB" in atom.fullname:
                        cb_coords.append(atom.get_coord())
                if residue.get_resname() == "GLY":
                    cb_coords.append([-10000000, -10000000, -10000000])
                residue_name_map[int(residue.get_id()[1])] = str(residue.get_resname())
                global_res_index[residue_counter] = int(residue.get_id()[1])
                residue_counter += 1
            break

    return (
        np.array(ca_coords),
        np.array(cb_coords),
        residue_name_map,
        global_res_index
        )


def read_pdbs_seed(active_residue_indices, seed_pdb_path):
    """
    Read a PDB file from the specified path, parse CA and CB coordinates, identify active site residues, 
    and return relevant data structures for further analysis.

    Parameters
    ----------
    active_residue_indices : list of int
        Indices of residues considered active in the structure.
    seed_pdb_path : str
        File path to the PDB file used as the seed structure.

    Returns
    -------
    numpy.ndarray
        Array of C-alpha coordinates from all residues in the seed PDB.
    numpy.ndarray
        Array of C-alpha coordinates from the active site residues.
    numpy.ndarray
        Array of C-beta coordinates from the active site residues.
    dict
        A dictionary mapping each residue ID (int) to its three-letter residue name (str).
    numpy.ndarray
        Array of active site residue IDs.
    dict
        Mapping of residue IDs to their corresponding array indices.

    Notes
    -----
    For 'GLY' residues, a placeholder coordinate [-10000000, -10000000, -10000000] 
    is inserted to represent the missing C-beta atom.
    """
    parser = PDBParser()
    parsed_seed_structure = parser.get_structure("complex", seed_pdb_path)

    seed_ca_coords = []
    active_site_ca_coords = []
    active_site_cb_coords = []
    seed_residue_names = {}
    active_site_seed_resids = []
    seed_resid_to_array_index = {}
    index_counter = 0

    for model in parsed_seed_structure:
        for chain in model:
            for residue in chain:
                seed_resid_to_array_index[int(residue.get_id()[1])] = index_counter
                index_counter += 1

                for atom in residue:
                    if "CA" in atom.fullname:
                        seed_ca_coords.append(atom.get_coord())

                seed_residue_names[int(residue.get_id()[1])] = str(residue.get_resname())

                if residue.get_id()[1] in active_residue_indices:
                    active_site_seed_resids.append(int(residue.get_id()[1]))

                    for atom in residue:
                        if "CA" in atom.fullname:
                            active_site_ca_coords.append(atom.get_coord())
                        if "CB" in atom.fullname:
                            active_site_cb_coords.append(atom.get_coord())

                    if residue.get_resname() == "GLY":
                        active_site_cb_coords.append([-10000000, -10000000, -10000000])

    return (
        np.array(seed_ca_coords),
        np.array(active_site_ca_coords),
        np.array(active_site_cb_coords),
        seed_residue_names,
        np.array(active_site_seed_resids),
        seed_resid_to_array_index
    )

def read_pdbs_seed_with_radius(active_residue_indices, seed_pdb_path, radius):
    """
    Read a PDB file from the specified path, parse CA and CB coordinates, identify active site residues, 
    and return relevant data structures for further analysis.

    Parameters
    ----------
    active_residue_indices : list of int
        Indices of residues considered active in the structure.
    seed_pdb_path : str
        File path to the PDB file used as the seed structure.

    Returns
    -------
    numpy.ndarray
        Array of C-alpha coordinates from all residues in the seed PDB.
    numpy.ndarray
        Array of C-alpha coordinates from the active site residues.
    numpy.ndarray
        Array of C-beta coordinates from the active site residues.
    dict
        A dictionary mapping each residue ID (int) to its three-letter residue name (str).
    numpy.ndarray
        Array of active site residue IDs.
    dict
        Mapping of residue IDs to their corresponding array indices.

    Notes
    -----
    For 'GLY' residues, a placeholder coordinate [-10000000, -10000000, -10000000] 
    is inserted to represent the missing C-beta atom.
    """
    parser = PDBParser()
    parsed_seed_structure = parser.get_structure("complex", seed_pdb_path)

    seed_ca_coords = []
    active_site_ca_coords = []
    active_site_cb_coords = []
    seed_residue_names = {}
    active_site_seed_resids = []
    seed_resid_to_array_index = {}
    index_counter = 0

    for model in parsed_seed_structure:
        for chain in model:
            for residue in chain:
                seed_resid_to_array_index[int(residue.get_id()[1])] = index_counter
                index_counter += 1

                for atom in residue:
                    if "CA" in atom.fullname:
                        seed_ca_coords.append(atom.get_coord())

                seed_residue_names[int(residue.get_id()[1])] = str(residue.get_resname())
                if residue.get_id()[1] in active_residue_indices:
                    active_site_seed_resids.append(int(residue.get_id()[1]))

                    for atom in residue:
                        if "CA" in atom.fullname:
                            active_site_ca_coords.append(atom.get_coord())
                        if "CB" in atom.fullname:
                            active_site_cb_coords.append(atom.get_coord())

                    if residue.get_resname() == "GLY":
                        active_site_cb_coords.append([-10000000, -10000000, -10000000])

    for model in parsed_seed_structure:
        for chain in model:
            for residue in chain:
                index_counter += 1

                for atom in residue:
                    if "CA" in atom.fullname:
                        dist = np.linalg.norm(active_site_ca_coords[0] - atom.get_coord())
                        break

                if dist < radius and dist > 0.02: 
                    active_site_seed_resids.append(int(residue.get_id()[1]))
                    seed_residue_names[int(residue.get_id()[1])] = "UNK"
                    for atom in residue:
                        if "CA" in atom.fullname:
                            active_site_ca_coords.append(atom.get_coord())
                        if "CB" in atom.fullname:
                            active_site_cb_coords.append(atom.get_coord())

                    if residue.get_resname() == "GLY":
                        active_site_cb_coords.append([-10000000, -10000000, -10000000])

    

    return (
        np.array(seed_ca_coords),
        np.array(active_site_ca_coords),
        np.array(active_site_cb_coords),
        seed_residue_names,
        np.array(active_site_seed_resids),
        seed_resid_to_array_index
    )

def print_protein(translation_vector, rotation, path, name):
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
    print("Saving the structure in: ", config.path_results+"/"+name+".pdb")
    io.save(config.path_results+"/"+name+".pdb")

def get_cavities(protein, res_active_site, volume_cutoff=5.0, removal_distance=2.4, probe_in=1.4, probe_out=4, surface="SES", step=0.6):
    """
    Retrieve cavity information from a given protein structure based on specified parameters.
    Args:
        protein (str): The file path to the target protein (in PDB format).
        res_active_site (list): A list of residues that define the active site of interest.
        volume_cutoff (float, optional): Minimum volume required for a cavity to be considered valid.
        removal_distance (float, optional): Distance threshold used to remove or merge cavities close to each other.
        probe_in (float, optional): Radius of the small probe for detecting cavity inner boundaries.
        probe_out (float, optional): Radius of the large probe for detecting cavity outer boundaries.
        surface (str, optional): The surface type to be considered in the detection (e.g., "SES").
        step (float, optional): Resolution step size for the internal volumetric calculation.
    Returns:
        list: A list of residue information that contributes to the detected cavities intersecting with the active site.
    """
    atomic = pyKVFinder.read_pdb(protein)
    volume_cutoff= 5.0
    removal_distance= 2.4
    probe_in= 1.4
    probe_out= 4
    surface="SES"
    step=0.6
    vertices = pyKVFinder.get_vertices(atomic, probe_out=probe_out, step=step)
    ncav, cavities = pyKVFinder.detect(atomic, vertices, step=step, probe_in=probe_in, probe_out=probe_out, removal_distance=removal_distance, volume_cutoff=volume_cutoff, surface=surface, nthreads=1)
    residues = pyKVFinder.constitutional(cavities, atomic, vertices, step=step, probe_in=probe_in, ignore_backbone=False, nthreads=1)
    
    inside=[]
    case_selected=[]
    for k, res in residues.items():
        selected=False
        for r in res:
            for aa in res_active_site:
                if str(aa) == r[0]:
                    selected =True
                    break
            if selected:
                break
        if selected:
            for r in res:
                if r[0] not in inside:
                    case_selected.append(r)
                    inside.append(r[0])
    return case_selected


def ActSeek_main(aa_des,  case_protein_filename, iterations, case_protein_name, seed_selected,seed_coords, cavity_coords,  aaCav, active, real_index_seed, cavity_coords_used, cavity_coords_cb_used, active_site):
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
        protein_coords, protein_coords_cb, aa, real_index= read_pdbs_case(case_protein_filename)
    except:
        return None, None, None, None

    
    #print(case_protein_filename)
    # Produces a vector with all the possible amino acid correspondences based on the aa_des dictionary
    pc = active_site.get_possible_correspondences(protein_coords, aa, real_index)
    threshold = config.threshold_others

    #print(pc)

    # Calculate the valid combinations by removing all the combinations of 3 amino acid mapping that are not in the appropiated distances with a threshold of 1 (that coud be an argument for the main function to be defined by the user)
    valid_combinations = active_site.get_all_possible_combinations(pc, protein_coords, aa, real_index, threshold,
                                                                   aa_des)

    #print(valid_combinations)

    if len(valid_combinations) > 0:
        distances, distances_arround, t_transformed, solution, translation_vector, rotation = ActSeekLib.actseek_search(
            cavity_coords, protein_coords, protein_coords_cb, seed_coords, valid_combinations,
            iterations, aa_des, aa, real_index, real_index_seed, aaCav, active, cavity_coords_used,cavity_coords_cb_used, config.aa_surrounding, config.threshold_combinations, 0)

        if np.sum(distances) / len(solution) < config.threshold and len(solution) >= 3 and distances_arround < config.aa_surrounding_threshold:
            try:
                rmsd, minrmsd, percentage,distancesal, indices= ActSeekLib.getGlobalDistance(t_transformed, seed_coords)

                if config.KVFinder == True:
                    try:
                        residues = []
                        for sol in solution:
                            residues.append(str(real_index[sol[0]]))
                        case_selected = get_cavities(case_protein_filename, residues) 

                        if seed_selected != None and case_selected != None:
                            kvdistance, kvmapping, kvpercentage = ActSeekLib.compare_cavity(case_selected, seed_selected, distancesal, indices)  
                            case_cavity=[str(res[0])+":"+res[2] for res in case_selected]  
                            
                    except:
                        raise
                        case_cavity=[]  

            except:
                traceback.print_exc()

            if len(config.testing) > 4:
                sol_write = open(f"{config.path_results}/{case_protein_name}.csv", "w")
                sol_write.write("Uniprot ID,Mapping,Average distance,Average distance AA arround, All distances,Structural local similarity, Structural RMSD, Percentage structural mapping,Cavity, Cavity distance, Cavity mapping (case:seed),Cavity mapping percentage\n")
                print_protein(translation_vector, rotation, case_protein_filename, case_protein_name)
            else:
                sol_write = open(f"{config.random_dir}/{case_protein_name}.txt", "w")

            strmapping = ""
            for aamap in solution:
                aa_protein = aa[real_index[aamap[0]]]
                aa_cavity = aaCav[active[aamap[1]]]
                strmapping = strmapping + str(real_index[aamap[0]]) + aa_protein + "-" + str(
                    active[aamap[1]]) + aa_cavity + ";"


            if config.KVFinder == True:                
                try:
                    sol_write.write(case_protein_name + "," + strmapping + "," + str(np.sum(distances) / len(solution)) + "," + str(distances_arround) + "," + ";".join(map(str, distances)) + ","+str(rmsd)+","+str(minrmsd)+","+str(percentage)+","+";".join(case_cavity)+","+str(kvdistance)+","+ ";".join(kvmapping)+","+ str(kvpercentage)+"\n")
                except:
                    sol_write.write(
                    case_protein_name + "," + strmapping + "," + str(np.sum(distances) / len(active)) + "," + str(distances_arround) + "," + ";".join(map(str, distances)) + ","+str(rmsd)+","+str(minrmsd)+","+str(percentage)+","+"Failed to calculate the cavity.\n")
            else:   
                sol_write.write(
                    case_protein_name + "," + strmapping + "," + str(np.sum(distances) / len(active)) + "," + str(distances_arround) + "," + ";".join(map(str, distances)) + ","+str(rmsd)+","+str(minrmsd)+","+str(percentage)+"\n")
            sol_write.close()
    else:        
        pass



def process_protein(case_protein_name,    
    seed_coords=None,
    cavity_coords=None,
    aaCav=None,
    active=None,
    real_index_seed=None,
    seed_selected=None,
    cavity_coords_cb_used = None,
    cavity_coords_used=None,
    active_site=None):
  
    """
    Processes a protein structure by downloading it from the AlphaFold database if needed, 
    running the ActSeek algorithm, and optionally removing the downloaded file when finished.
    Args:
        case_protein_name (str): The protein identifier used to construct the filename 
            for AlphaFold-based retrieval.
        seed_coords (any, optional): Seed coordinates used by the ActSeek algorithm (default: None).
        cavity_coords (any, optional): Cavity coordinates for the ActSeek algorithm (default: None).
        aaCav (any, optional): Amino acid cavity data (default: None).
        active (any, optional): Indicator for the active state (default: None).
        real_index_seed (any, optional): Actual index for the seed (default: None).
        seed_selected (any, optional): Selection criteria for the seed (default: None).
        cavity_coords_cb_used (any, optional): Specific cavity coordinates used in calculations (default: None).
        cavity_coords_used (any, optional): Cavity coordinates used in computations (default: None).
        active_site (any, optional): Information about the active site (default: None).
    Raises:
        requests.exceptions.RequestException: If an error occurs during the retrieval 
            of the protein structure from the AlphaFold database.
        OSError: If there's an issue removing the downloaded protein file.
    Note:
        This function depends on global configuration (e.g., config.alphafold_proteins_path, 
        config.iterations, config.delete_protein_files) as well as external libraries 
        (like requests, wget, and the internal ActSeek_main function).
    """
    case_protein_name = case_protein_name

    try:
        if config.custom == False:
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

        else:
            case_protein_filename = case_protein_name+".pdb"
            case_protein_filepath = f"{config.alphafold_proteins_path}/{case_protein_filename}" 

        # Main function of the algorithm       
        ActSeek_main(config.aa_grouping, case_protein_filepath, config.iterations, case_protein_name, seed_selected,seed_coords, cavity_coords,aaCav, active, real_index_seed, cavity_coords_used, cavity_coords_cb_used,active_site)

        # Removes the protein structure once the algorithm has finished
        if config.delete_protein_files == True:
            try:
                os.remove(case_protein_filepath)
            except:
                traceback.print_exc()

                pass
    except:
        pass


def read_config(file_path):
    with open(file_path, 'r') as file:
        config = json.load(file)
    return config




def parse_args():
    """
    Parse the command line arguments and return them as a dictionary.
    """
    try:
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
        parser.add_argument("-kv", "--KVFinder",dest="KVFinder", action="store_true", default=config['KVFinder'],
                        help="Uses KVFinder to compare the cavity where the active sides are with the cavity in the seed structure.")
        parser.add_argument("-c", "--custom",dest="custom", action="store_true", default=config['custom'],
                        help="Using a custom structure database not connected with Alphafold for the search.")
        parser.add_argument("-rd", "--radius", type=float, default=config['radius'],
                            help="Threshold of the average distance of the mapped amino acids.")
        
        
        args = parser.parse_args()
        return vars(args)
    except:
        print("Config file not found.")
        return None

def main():
    """
    Main function to set up configuration, reads and process the seed protein, create necessary directories,
    and process proteins using multiprocessing.
    """
   
    global config
    config = parse_args()   

    if config == None:
        return 
    
    if os.path.isdir(config['path_results']) == False:
        os.mkdir(config['path_results'])

    if os.path.isdir(config['alphafold_proteins_path']) == False:
        os.mkdir(config['alphafold_proteins_path'])
        
    with tempfile.TemporaryDirectory(dir=config['path_results']) as tmpdir:
        config['random_dir'] = tmpdir
        

        config = argparse.Namespace(**config)

        # Set random seed to a fix value, to ensure repeatability
        random.seed(config.random_seed)


        seed_protein = config.seed_protein_file
        aa_active = [int(x) for x in config.active_site.split(",")]
        #print(aa_active)

        if len(aa_active) < 3: 
            try:
                radius = config.radius
                seed_coords, cavity_coords, cavity_coords_cb, active_amino_acids_dict, active_amino_acid_index, real_index_seed_dict = read_pdbs_seed_with_radius(aa_active,
                    seed_protein, radius)                 
            except:
                print("Seed pdb file not found or the file is not readable.")
                return
            amino_acids_used_for_search = np.array([0,1,2])
        else:
            try:
                seed_coords, cavity_coords, cavity_coords_cb, active_amino_acids_dict, active_amino_acid_index, real_index_seed_dict = read_pdbs_seed(aa_active,
                    seed_protein)
            except:
                print("Seed pdb file not found or the file is not readable.")
                return
        
            amino_acids_used_for_search = np.array([int(x) for x in config.selected_active.split(",")])

        active_used_for_search = np.array([active_amino_acid_index[x] for x in amino_acids_used_for_search])
        active_coords_used_for_search = np.array([cavity_coords[x] for x in amino_acids_used_for_search])
        active_coords_cb_used_for_search = np.array([cavity_coords_cb[x] for x in amino_acids_used_for_search])

        #print(active_used_for_search)
        # Creates an object of the class Active_side where all the information is added
        active_site = ActSeekLib.Active_site(active_used_for_search, active_amino_acids_dict, active_coords_used_for_search, config.aa_grouping, config.radius)
        active_site.get_AA_groups(False)
        

        
        '''residues = active_site.amino_acid_list
        print(residues[0].index)
        print(residues[0].real_index)        
        print(residues[0].amino_acid_class)
        print(residues[0].amino_acid_id)
        print(residues[0].distances)


        print(residues[1].index)
        print(residues[1].real_index)        
        print(residues[1].amino_acid_class)
        print(residues[1].amino_acid_id)
        print(residues[1].distances)


        print(residues[2].index)
        print(residues[2].real_index)        
        print(residues[2].amino_acid_class)
        print(residues[2].amino_acid_id)
        print(residues[2].distances)'''



        if config.KVFinder == True:
            seed_selected = get_cavities(seed_protein, active_amino_acid_index) 
            
        else:
            seed_selected=None

       

        processProteinWithData = partial(
            process_protein,
            seed_coords=seed_coords,
            cavity_coords=cavity_coords,
            aaCav=active_amino_acids_dict,
            active=active_amino_acid_index,
            real_index_seed=real_index_seed_dict,
            seed_selected=seed_selected,
            cavity_coords_cb_used = active_coords_cb_used_for_search,
            cavity_coords_used=active_coords_used_for_search,
            active_site=active_site
        )

       
        if len(config.testing) > 4:
            case_protein = config.testing
            processProteinWithData(case_protein)
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
                    futures = [executor.submit(processProteinWithData, case) for case in cases]
                    for future  in concurrent.futures.as_completed(futures):
                        pbar.update(1)

            files = os.listdir(config.random_dir)
            results = open(config.path_results+"/results.csv","w")
            results.write("Uniprot ID,Mapping,Average distance,Average distance AA arround, All distances,Structural local similarity, Structural RMSD, Percentage structural mapping,Cavity, Cavity distance, Cavity mapping (case:seed),Cavity mapping percentage\n")
            for file in files:
                f = open(config.random_dir+"/"+file, "r")
                for line in f:
                    results.write(line)
                f.close()
                os.remove(config.random_dir+"/"+file)
            results.close()


if __name__ == '__main__':
    main()
