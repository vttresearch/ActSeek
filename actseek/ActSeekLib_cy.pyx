import itertools
import numpy as np
cimport numpy as cnp
import random
import warnings
from scipy.spatial import KDTree
import traceback
from actseek.grid import detect, get_vertices,spatial, constitutional, _get_sincos, _get_dimensions


warnings.filterwarnings("ignore")

config = None



cdef class Amino_acid:
    """
    This class keeps information about the amino acids in the active site. The information needed is:
    index: the index of this amino acid in the coordinates vector
    real_index: the corresponding index in the pdb file.
    amino_acid_id: amino acid identity (Ser, Pro..)
    amino_acid_class: it may be the same as the amino_acid_id or it could be different. The classes are defined in the dictionary "aa_des"
    coords: a numpy array with 3 float numbers corresponding to the amino acid 3D coordinates
    distances: a dictionary with the distance of this amino acid to all the amino acids of specific class.The key is the class and the item a vector with distances.
    """
    cdef public int index
    cdef public int real_index
    cdef public str amino_acid_class
    cdef public str amino_acid_id
    cdef public cnp.ndarray coords
    cdef public dict distances

    def __init__(self, int index, int real_index, str amino_acid_id, str amino_acid_class, cnp.ndarray coords):
        self.index = index
        self.real_index = real_index
        self.amino_acid_class = amino_acid_class
        self.amino_acid_id = amino_acid_id
        self.coords = coords
        self.distances = {}

    def get_coords(self):
        """
        Get the amino acid 3D coordinate
        """
        return self.coords

    def get_class(self):
        """
        Get the amino acid class. It can be the amino acid identity or the class of amino acid (the classes are defined in the variable aa_des)
        """
        return self.amino_acid_class

    def set_distances(self, coord, amino_acid_class):
        """
        populates the distances dictionary with the distance from this amino acid to other amino acids of different classes.
        """
        if amino_acid_class in self.distances:
            dist = self.distances[amino_acid_class]
            dist.append(np.linalg.norm(self.coords - coord))
            self.distances[amino_acid_class] = dist
        else:
            dist = [np.linalg.norm(self.coords - coord)]
            # print(self.aa_class,aa_class, dist)
            self.distances[amino_acid_class] = dist

    def get_distance(self, amino_acid_class):
        """
        Get the distances from this amino acid to other amino acids of specific class.
        """
        if amino_acid_class in self.distances:
            return self.distances[amino_acid_class]
        else:
            return None

    def get_by_index(self):
        """
        Returns the amino acid index
        """
        return self.index


cdef class Active_site:
    """
    Class that keeps information about the amino acids that are involved in the search (usually active or binding sites)
    In the init function the list 'amino_acid_list' is populated with objects of class Amino_acids corresponding to the ones used for the search. 
    It also calculate the distances between the amino acids used in the search.    
    """
    cdef public cnp.ndarray search_list
    cdef public dict amino_acid_dict
    cdef public list amino_acid_list
    cdef public dict amino_acid_groups
    cdef public cnp.ndarray coords
    
    def __init__(self, search_list, amino_acid_dict, coords, amino_acid_groups):
        self.search_list = np.array(search_list)
        self.amino_acid_dict = amino_acid_dict
        self.amino_acid_list = []
        self.coords = coords
        self.amino_acid_groups = amino_acid_groups
        self.add_amino_acids()
        self.set_self_distance()

    def add_amino_acids(self):
        cdef int amino_acid_index
        cdef str amino_acid_id
        cdef object amino_acid
        cdef cnp.ndarray index
        for amino_acid_index, amino_acid_id in self.amino_acid_dict.items():
            index, = np.where(self.search_list == amino_acid_index)
            if len(index)>0:
                try:
                    amino_acid = Amino_acid(index[0], amino_acid_index, amino_acid_id, self.amino_acid_groups[amino_acid_id], self.coords[index[0]])
                    self.amino_acid_list.append(amino_acid)
                except:
                    j=0


    def get_Number_AA(self):        
        return len(self.amino_acid_list)

    def get_AA(self):
        return self.amino_acid_list

    def get_AA_groups(self, bint unique):
        cdef list groups
        cdef object amino_acid
        groups = []

        for amino_acid in self.amino_acid_list:
            groups.append(amino_acid.get_class())
        if unique:
            return list(set(groups))
        else:
            return groups

    def get_AAs_of_class(self, str amino_acid_class):
        cdef list class_list
        cdef object amino_acid
        cdef str group
        class_list = []
        for amino_acid in self.amino_acid_list:
            group = amino_acid.get_class()
            if group == amino_acid_class:
                class_list.append(amino_acid)

    def set_self_distance(self):
        cdef object amino_acid
        cdef object amino_acid2
        for amino_acid in self.amino_acid_list:
            for amino_acid2 in self.amino_acid_list:
                if amino_acid != amino_acid2:
                    amino_acid.set_distances(amino_acid2.coords, amino_acid2.get_class())

    def get_distances(self, str amino_acid_class, str amino_acid_class2):
        cdef list distances_list
        cdef object amino_acid
        cdef list distances
        distances_list = []
        for amino_acid in self.amino_acid_list:
            if amino_acid.get_class() == amino_acid_class:
                distances = amino_acid.get_distance(amino_acid_class2)
                for d in distances:
                    distances_list.append(d)
        return distances_list

    def get_aa_by_index(self, int index):
        cdef object amino_acid
        for amino_acid in self.amino_acid_list:
            if amino_acid.get_by_index() == index:
                return amino_acid

    def get_possible_correspondences(self, cnp.ndarray case_protein, dict case_protein_amino_acids, dict real_index):
        """
        It finds all the amino acids from the case protein that belong of the same classes represented in the list of amino acids used for the search. 
        """
        cdef list correspondences
        cdef int amino_acid1
        cdef object amino_acid2
        correspondences = []
        for amino_acid1 in range(0, len(case_protein)):
            for amino_acid2 in self.amino_acid_list:
                try:
                    if self.amino_acid_groups[case_protein_amino_acids[real_index[amino_acid1]]] == amino_acid2.get_class():
                        correspondences.append([amino_acid1, amino_acid2.get_by_index()])
                except:
                    pass


        return correspondences

    def check_distances(self, float distance1, list distances, float threshold):
        cdef float distance

        if len(distances) == 0:
            return False
        for distance in distances:
            try:
                if abs(distance1 - distance) <= threshold:
                    return True
            except:
                return False
        return False

    def calculate_distance_in_case_protein(self, int index1, int index2, cnp.ndarray protein_coords):
        return np.linalg.norm(protein_coords[index1] - protein_coords[index2])



    def get_all_possible_combinations(self, list correspondences, cnp.ndarray protein_coords, dict case_protein_amino_acid_dict, dict real_index, float threshold, dict amino_acid_groups):
        """
            Using the correspondences calculated in "get_possible_correspondences", it creates combinations of 3. 
            These combinations should have repeated indexes and the distances of the amino acis should be in the same distances with respect each other as the three amino acids used for the search.
        """
        cdef object combinations
        cdef list valid_combinations
        cdef tuple combo
        cdef float dist1
        cdef float dist2
        cdef float dist3
        cdef list distances1
        cdef list distances2
        cdef list distances3

        combinations = itertools.combinations(correspondences, 3)
        valid_combinations = []

        for combo in combinations:
            distinct_classes = set(item[1] for item in combo)
            if len(distinct_classes) < 3:
                continue

            dist1 = self.calculate_distance_in_case_protein(combo[0][0], combo[1][0], protein_coords)
            if dist1 is None:
                continue
            dist1 = self.calculate_distance_in_case_protein(combo[0][0], combo[1][0], protein_coords)
            if dist1 == None:
                continue

            distances1 = self.get_distances(amino_acid_groups[case_protein_amino_acid_dict[real_index[combo[0][0]]]], amino_acid_groups[case_protein_amino_acid_dict[real_index[combo[1][0]]]])
            if not self.check_distances(dist1, distances1, threshold):
                continue

            dist2 = self.calculate_distance_in_case_protein(combo[1][0], combo[2][0], protein_coords)
            if dist2 == None:
                continue

            distances2 = self.get_distances(amino_acid_groups[case_protein_amino_acid_dict[real_index[combo[1][0]]]], amino_acid_groups[case_protein_amino_acid_dict[real_index[combo[2][0]]]])
            if not self.check_distances(dist2, distances2, threshold):
                continue

            dist3 = self.calculate_distance_in_case_protein(combo[2][0], combo[0][0], protein_coords)
            if dist3 == None:
                continue

            distances3 = self.get_distances(amino_acid_groups[case_protein_amino_acid_dict[real_index[combo[2][0]]]], amino_acid_groups[case_protein_amino_acid_dict[real_index[combo[0][0]]]])
            if not self.check_distances(dist3, distances3, threshold):
                continue

            valid_combinations.append(combo)

        return valid_combinations



cpdef euclidean_transformation_fit(cnp.ndarray search_alphac_coords, cnp.ndarray search_betac_coords, cnp.ndarray case_protein_alphac_coords, cnp.ndarray case_protein_betac_coords, list combinations, int index):
    cdef tuple selected_combination
    cdef list search_coords
    cdef list case_protein_coords
    cdef list pair
    cdef cnp.ndarray rotation
    cdef cnp.ndarray translation

    if index == -1:
        selected_combination = tuple(random.sample(combinations, 1))
        selected_combination = selected_combination[0]
    else:
        selected_combination = combinations[index]
    search_coords = []
    case_protein_coords = []
    # take coordinates of more atoms
    for pair in selected_combination:
        search_coords.append(search_alphac_coords[pair[1]])
        case_protein_coords.append(case_protein_alphac_coords[pair[0]])
        if search_betac_coords[pair[1]][0] != -10000000 and case_protein_betac_coords[pair[0]][0] != -10000000:
            search_coords.append(search_betac_coords[pair[1]])
            case_protein_coords.append(case_protein_betac_coords[pair[0]])

    translation_vector, rotation = estimate_alignment(case_protein_coords, search_coords)

    return translation_vector, rotation, list(selected_combination)


# Euclidean transformation fitting
cpdef estimate_alignment(list points1,list points2):

    # Center the point sets
    cdef cnp.ndarray centroid1
    cdef cnp.ndarray centroid2
    cdef cnp.ndarray centered_points1
    cdef cnp.ndarray centered_points2

    cdef cnp.ndarray covariance_matrix
    cdef cnp.ndarray u
    cdef cnp.ndarray vh
    cdef cnp.ndarray _
    cdef cnp.ndarray rotation
    cdef cnp.ndarray translation

    try:
        centroid1 = np.mean(points1, axis=0)
        centroid2 = np.mean(points2, axis=0)
        centered_points1 = points1 - centroid1
        centered_points2 = points2 - centroid2

        # Compute the covariance matrix
        covariance_matrix = centered_points2.T @ centered_points1

        # Perform singular value decomposition (SVD)
        u, _, vh = np.linalg.svd(covariance_matrix)

        # Compute the rotation matrix
        rotation = vh.T @ u.T
        
        if np.linalg.det(rotation) < 0:
            vh.T[2] = -vh.T[2]
            rotation = np.transpose(np.dot(np.transpose(vh.T), np.transpose(u)))
            
        # Compute the translation vector
        translation = centroid2 - centroid1 @ rotation

        return translation, rotation
    except:
        return None,None




cpdef calculate_distances(cnp.ndarray transformed_case_protein_coords, cnp.ndarray search_alphac_coords, threshold):
    """
    Calculate distances between transformed case protein coordinates and search alpha carbon coordinates,
    ensuring each search alpha carbon coordinate has a unique nearest neighbor.

    Parameters:
    transformed_case_protein_coords (cnp.ndarray): Transformed case protein coordinates.
    search_alphac_coords (cnp.ndarray): Search alpha carbon coordinates.
    threshold (float): Distance threshold.

    Returns:
    tuple: Three elements -
        final_distances (np.ndarray): Array of distances below the threshold.
        final_indices (np.ndarray): Array of index pairs (transformed_case_protein_coords index, search_alphac_coords index).        
    """
    cdef cnp.ndarray distances
    cdef cnp.ndarray indices
    cdef list final_distances
    cdef list final_indices

    tree = KDTree(transformed_case_protein_coords)
    distances, indices = tree.query(search_alphac_coords)
    final_distances=[]
    final_indices=[]
    for i in range(len(distances)):
        if distances[i] < threshold:
            final_distances.append(distances[i])
            final_indices.append([indices[i], i])

    return np.array(final_distances), np.array(final_indices)


cpdef calculate_distance_single(int index1, int index2, cnp.ndarray case_protein_coords):
    return np.linalg.norm(case_protein_coords[index1] - case_protein_coords[index2])




def calculate_final_distance(search_alphac_coords_all, case_protein_alphac_coords, case_protein_betac_coords, valid_combinations, index, amino_acid_groups, case_protein_dict, real_index, search_amino_acids_dict, search_indexes, search_alphac_coords, search_betac_coords, threshold):
    """
    Calculate the final distance between the transformed case protein coordinates and the search coordinates. The variable search_alphac_coords_all is different from search_alphac_coords as it contains all the amino acids defined by the user,
    while search_alphac_coords contains only the 3 selected amino acids used for the actual search.

    Parameters:
    -----------
    search_alphac_coords_all : array-like
        Contains all the amino acids defined by the user.
    case_protein_alphac_coords : array-like
        Alpha carbon coordinates of the case protein.
    case_protein_betac_coords : array-like
        Beta carbon coordinates of the case protein.
    valid_combinations : array-like
        Valid combinations of indices for the transformation.
    index : int
        Index of the current combination.
    amino_acid_groups : dict
        Dictionary mapping amino acids to their groups.
    case_protein_dict : dict
        Dictionary mapping indices to amino acids in the case protein.
    real_index : array-like
        Real indices of the case protein.
    search_amino_acids_dict : dict
        Dictionary mapping indices to amino acids in the search set.
    search_indexes : array-like
        Indices of the search amino acids.
    search_alphac_coords : array-like
        Alpha carbon coordinates of the selected search amino acids.
    search_betac_coords : array-like
        Beta carbon coordinates of the selected search amino acids.
    threshold : float
        Distance threshold for valid mappings.

    Returns:
    --------
    mapping : np.ndarray
        Array of mapped indices between the case protein and search coordinates.
    distance_sum : float
        Sum of the distances for the valid mappings.
    distances : list
        List of distances for the valid mappings.
    t_transformed : array-like
        Transformed alpha carbon coordinates of the case protein.
    translation_vector : array-like
        Translation vector used for the transformation.
    rotation : array-like
        Rotation matrix used for the transformation.
    """
    

    translation_vector, rotation, mapsel = euclidean_transformation_fit(search_alphac_coords, search_betac_coords, case_protein_alphac_coords, case_protein_betac_coords,
                                                                        valid_combinations, index)

    if rotation is None:
        return None, None, None, None, None, None

    t_transformed = case_protein_alphac_coords @ rotation + translation_vector

    dist, indices = calculate_distances(t_transformed, search_alphac_coords_all, threshold)

    distance_sum = 0
    mapping = []
    distances = []
    amino_acid_groups_case = [amino_acid_groups[case_protein_dict[real_index[ind[0]]]] for ind in indices]
    amino_acid_groups_search = [amino_acid_groups[search_amino_acids_dict[search_indexes[e]]] for e in range(len(search_alphac_coords_all))]

    for e in range(len(search_alphac_coords_all)):
        valid_indices = [i for i, ind in enumerate(indices) if ind[1] == e and amino_acid_groups_case[i] == amino_acid_groups_search[e]]
        if valid_indices:
            minindex = min(valid_indices, key=lambda i: dist[i])
            distances.append(dist[minindex])
            distance_sum += dist[minindex]
            mapping.append(indices[minindex])

    if not mapping:
        distance_sum = 100000

    return np.array(mapping), distance_sum, distances, t_transformed, translation_vector, rotation




cpdef get_distance_around(cnp.ndarray case_protein_coords, cnp.ndarray search_protein_coords, cnp.ndarray mapping, dict search_real_index, cnp.ndarray search_indexes, int number_of_surrounding, int printing):
    """
    Calculate the sum and mean of distances between "number_of_surrounding" number of corresponding amino acids in two protein structures surrounding the amino acids selected for the search.

    Parameters
    ----------
    case_protein_coords : cnp.ndarray
        Coordinates of the case protein.
    search_protein_coords : cnp.ndarray
        Coordinates of the search protein.
    mapping : cnp.ndarray
        Mapping of amino acid indices between the case and search proteins.
    search_real_index : dict
        Dictionary mapping search protein indices to real indices.
    search_indexes : cnp.ndarray
        Array of search protein indices.
    number_of_surrounding : int
        Number of surrounding amino acids to consider for distance calculation.
    printing : int
        Flag to enable or disable printing of intermediate results (1 to enable, 0 to disable).

    Returns
    -------
    tuple
        A tuple containing the sum of distances and the mean distance.
    """

    cdef list distances = []
    cdef int i, case_protein_index, search_index
    cdef float distance, mean_distances

    for amino_acids_mapped in mapping:
        case_protein_index = amino_acids_mapped[0]
        search_index = search_real_index[search_indexes[amino_acids_mapped[1]]]
        for i in range(1, number_of_surrounding + 1):
            if case_protein_index - i >= 0 and search_index - i >= 0:
                distance = np.linalg.norm(case_protein_coords[case_protein_index - i] - search_protein_coords[search_index - i])
                distances.append(distance)
                if printing == 1:
                    print("minus", i, case_protein_index - i, search_index - i, case_protein_coords[case_protein_index - i], search_protein_coords[search_index - i], distance)

            if case_protein_index + i < len(case_protein_coords) and search_index + i < len(search_protein_coords):
                distance = np.linalg.norm(case_protein_coords[case_protein_index + i] - search_protein_coords[search_index + i])
                distances.append(distance)
                if printing == 1:
                    print("plus", i, case_protein_index + i, search_index + i, case_protein_coords[case_protein_index + i], search_protein_coords[search_index + i], distance)

    mean_distances = np.mean(distances) if distances else 0.0
    return np.sum(distances), mean_distances




def ransac_protein(search_alphac_coords_all, case_protein_alphac_coords, case_protein_betac_coords, search_protein_coords, valid_combinations, iterations, amino_acids_groups, case_protein_dict,
                real_index, search_real_index, search_amino_acids_dict, search_indexes, search_alphac_coords, search_betac_coords, number_of_surrounding, threshold, printing):
    """
    Performs RANSAC (Random Sample Consensus) algorithm to find the best alignment of a case protein to a search protein.

    Parameters:
        search_alphac_coords_all (list): List of alpha carbon coordinates for the search protein.
        case_protein_alphac_coords (list): List of alpha carbon coordinates for the case protein.
        case_protein_betac_coords (list): List of beta carbon coordinates for the case protein.
        search_protein_coords (list): List of coordinates for the search protein.
        valid_combinations (list): List of valid combinations of amino acids.
        iterations (int): Number of iterations for the RANSAC algorithm.
        amino_acids_groups (dict): Dictionary of amino acid groups.
        case_protein_dict (dict): Dictionary of case protein information.
        real_index (list): List of real indices for the case protein.
        search_real_index (list): List of real indices for the search protein.
        search_amino_acids_dict (dict): Dictionary of search protein amino acids.
        search_indexes (list): List of search protein indices.
        search_alphac_coords (list): List of alpha carbon coordinates for the search protein.
        search_betac_coords (list): List of beta carbon coordinates for the search protein.
        number_of_surrounding (int): Number of surrounding residues to consider.
        threshold (float): Distance threshold for considering a match.
        printing (bool): Flag to enable or disable printing of debug information.

    Returns:
        tuple: A tuple containing:
            - distances_selected (list): List of selected distances.
            - final_distances_arround (float): Final mean distance around the selected mapping.
            - case_protein_transformed_final (list): Transformed coordinates of the case protein.
            - solution (list): List of mappings between case and search proteins.
            - translation_matrix (list): Translation vector for the transformation.
            - rotation_matrix (list): Rotation matrix for the transformation.
    """
    max_distance = float('inf')
    solution = []
    distances_selected = []
    translation_matrix = None
    rotation_matrix = None
    case_protein_transformed_final = None
    final_distances_arround = None

    num_combinations = len(valid_combinations)
    num_iterations = min(num_combinations, iterations)

    for i in range(num_iterations):
        try:
            j = -1 if num_combinations > iterations else i

            mapping, distances_sum, distances, case_protein_transformed, translation_vector, rotation = calculate_final_distance(
                search_alphac_coords_all, case_protein_alphac_coords, case_protein_betac_coords, valid_combinations, j, amino_acids_groups, case_protein_dict, real_index, search_amino_acids_dict, search_indexes, search_alphac_coords, search_betac_coords, threshold)

            if rotation is None:
                continue

            mpn = [mp[0] for mp in mapping]

            distances_arround_sum, distances_arround_mean = get_distance_around(
                case_protein_transformed, search_protein_coords, mapping, search_real_index, search_indexes, number_of_surrounding, printing)

            total_distance_sum = distances_sum + distances_arround_sum
            avg_distance = total_distance_sum / len(mapping)

            if avg_distance <= max_distance:
                max_distance = avg_distance
                final_distances_arround = distances_arround_mean
                solution = mapping
                distances_selected = distances
                translation_matrix = translation_vector
                rotation_matrix = rotation
                case_protein_transformed_final = case_protein_transformed

            if max_distance == 0 and len(set(mpn)) == len(search_alphac_coords_all):
                return distances_selected, distances_arround_mean, case_protein_transformed_final, solution, translation_matrix, rotation_matrix

        except Exception:
            traceback.print_exc()

    return distances_selected, final_distances_arround, case_protein_transformed_final, solution, translation_matrix, rotation_matrix
def ransac_protein(search_alphac_coords_all, case_protein_alphac_coords, case_protein_betac_coords, search_protein_coords, valid_combinations, iterations, amino_acids_groups, case_protein_dict,
                real_index, search_real_index, search_amino_acids_dict, search_indexes, search_alphac_coords,search_betac_coords, number_of_surrounding, threshold, printing):

    max_distance = 100000000
    solution = []
    distances_selected = []
    translation_matrix = None
    rotation_matrix = None
    case_protein_transformed_final = None
    distances_arround_mean= None
    final_distances_arround=None

    # If the number of valid combinations is less than the number of iterations, the algorithm test each combination one by one, if not, it takes the combination randomly
    for i in range(len(valid_combinations)):
        try:
            if len(valid_combinations)> iterations:
                j = -1
                if i > iterations:
                    break
            else:
                j=i

            mapping, distances_sum, distances, case_protein_transformed, translation_vector, rotation= calculate_final_distance(search_alphac_coords_all, case_protein_alphac_coords, case_protein_betac_coords, valid_combinations, j, amino_acids_groups, case_protein_dict, real_index, search_amino_acids_dict, search_indexes, search_alphac_coords,search_betac_coords, threshold)
            if type(rotation) == type(None):
                continue

            # Select the mappings with less distances. There are few cases where the protein has 2 or more possible mappings and the one with less distance is not the right one

            mpn = [mp[0] for mp in mapping]

            distances_arround_sum, distances_arround_mean = get_distance_around(case_protein_transformed, search_protein_coords,
                                                                                mapping, search_real_index, search_indexes,number_of_surrounding, printing)

            distances_sum = distances_sum + distances_arround_sum
            if (distances_sum/len(mapping))/len(mapping) <= max_distance:
                max_distance = (distances_sum/len(mapping))/len(mapping)

                final_distances_arround = distances_arround_mean
                solution = mapping
                distances_selected = distances
                translation_matrix = translation_vector
                rotation_matrix = rotation
                case_protein_transformed_final = case_protein_transformed

            if max_distance == 0 and len(set(mpn)) == len(search_alphac_coords_all):
                return distances_selected, distances_arround_mean, case_protein_transformed_final, solution, translation_matrix, rotation_matrix

        except Exception:
            traceback.print_exc()


    return  distances_selected, final_distances_arround, case_protein_transformed_final, solution, translation_matrix, rotation_matrix



def get_atomics(protein_coords,protein_coords_cb, protein_amino_acid_dict, real_index):
    atomics = []
    atoms_size ={"GLY": "2", "ALA": "2", "PRO": "2.2", "ARG": "5", "HIS": "4.5", "LYS": "4.75", "ASP": "4",
                  "GLU": "4.2",
                  "SER": "2.5", "THR": "2.4", "ASN": "3.5", "GLN": "4.2", "CYS": "2.5", "VAL": "2", "ILE": "2.5",
                  "LEU": "2.5",
                  "MET": "3", "PHE": "4", "TYR": "4", "TRP": "5"}
    for i in range(len(protein_coords)):
        coords = protein_coords[i]        
        coords2 = protein_coords_cb[i]
        res =  protein_amino_acid_dict[real_index[i]]
        atomics.append( [real_index[i],"A",res, "C", coords[0], coords[1], coords[2], 3])
        if coords2[0]> -1000: 
            atomics.append( [real_index[i],"A", res, "C", coords2[0], coords2[1], coords2[2], atoms_size[res]])
    
    return np.asarray(atomics)


def get_cavity(atomic):

    ligand_cutoff: Union[float, in32t] = 6
    volume_cutoff: Union[float, int32] = 2.0
    removal_distance: Union[float, int32] = 3
    probe_in: Union[float, int32] = 1.4
    probe_out: Union[float, int32] = 5

    step=0.6
    vertices = get_vertices(atomic,step,1)


    dim = _get_dimensions(vertices, step)
    nx = dim[0]
    ny = dim[1]
    nz = dim[2]
    rotation = _get_sincos(vertices)


    

    ncav, cavities = detect(
        atomic,
        vertices,
        step,
        probe_in,
        probe_out,
        removal_distance,
        volume_cutoff,
        None,
        ligand_cutoff,
        True,
        "SES",
        1,
        False,
    )
    if ncav > 0:
        residues = constitutional(cavities,atomic, vertices, ignore_backbone=False)
        return residues
    else:
        return []




cpdef find_nearest_neighbors(cnp.ndarray vector_source, cnp.ndarray vector_target, float threshold):
    """
    Find the nearest neighbors between two sets of vectors within a given threshold.

    Parameters
    ----------
    vector_source : cnp.ndarray
        The source vector set.
    vector_target : cnp.ndarray
        The target vector set.
    threshold : float
        The distance threshold for considering neighbors.

    Returns
    -------
    tuple
        A tuple containing:
        - distances (list): List of distances between aligned vectors.
        - alignment (list): List of index pairs (source index, target index) for aligned vectors.
    """
    m = len(vector_source)
    n = len(vector_target)
    
   
    dist_matrix = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            dist = np.linalg.norm(np.array(vector_source[i]) - np.array(vector_target[j]))
            dist_matrix[i][j] = dist
    
   
    score_matrix = (dist_matrix <= threshold).astype(int)
    
   
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
 
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            current_score = score_matrix[i-1][j-1]
            diagonal = dp[i-1][j-1] + current_score
            up = dp[i-1][j]
            left = dp[i][j-1]
            dp[i][j] = max(diagonal, up, left)
    
    i, j = m, n
    alignment = []
    distances = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + score_matrix[i-1][j-1]:
            if score_matrix[i-1][j-1] == 1:
                alignment.append((i-1, j-1))
                distances.append(dist_matrix[i-1][j-1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j]:
            i -= 1
        else:
            j -= 1    
    
    alignment.reverse()
    return distances, alignment

cpdef score(cnp.ndarray vector, int size):
    """
    Scores the given vector.

    Parameters
    ----------
    vector : cnp.ndarray
        The input vector to be scored.
    size : int
        The size of the input vector.

    Returns
    -------
    float
        The score of the input vector.
    """

    cdef float consecutive_count
    cdef float score
    max_consecutive_chunk = 0
    chunk_count = 0
    current_chunk = 1

    for i in range(1, len(vector)):
        if vector[i] == vector[i - 1] + 1:
            current_chunk += 1
        if current_chunk > size:
            chunk_count += 1
            current_chunk = 1

    # Check the last chunk
    if current_chunk > size:
        chunk_count += 1
    score = (chunk_count*size) / len(vector)
    return score



cpdef float score_vector(list vector):
    """
    Calculate a weighted score for the given vector.

    Parameters
    ----------
    vector : list
        The input vector for which the score is to be calculated.

    Returns
    -------
    float
        The weighted score of the input vector.

    Notes
    -----
    The function calculates the weighted score by iterating over a predefined
    set of weights [5, 10, 15, 20]. For each weight, it calls the `score` function
    with the input vector and the current weight and adds it to the total score. The 
    returned value is the weighted average of the scores.
    """
    cdef float total_score = 0.0
    cdef int total_weight = 0
    cdef int i
    target_indices=[]
    for val in vector:
        target_indices.append(val[1])
    target_indices = np.array(target_indices)
    for i in [5, 10, 15, 20]:
        total_score += score(target_indices, i) * i
        total_weight += i

    return total_score / total_weight

cpdef getGlobalDistance(cnp.ndarray coords1, cnp.ndarray coords2):
    """
    This function calculates the nearest neighbors between two sets of coordinates, computes a score based on the indices of the nearest neighbors, and returns the score, the root mean square deviation (RMSD) normalized by coverage, and the coverage.

    Returns:
        tuple: A tuple containing the score, normalized RMSD, and coverage. If an exception occurs, it returns (40.0, 100.0, 100.0).

    Raises:
        Exception: If an error occurs during the computation, the exception is caught, and a traceback is printed.
    """
    cdef list distances
    cdef list indices
    cdef float rmsd
    cdef float score

    try:
        distances, indices= find_nearest_neighbors(coords1, coords2,2)        
        score = score_vector(indices)


        if len(indices) == 0:
            return 100.0, 100.0, 100.0,[],[]

        rmsd = np.mean(distances)
        coverage = len(indices) / len(coords1)

        return score, rmsd / coverage, coverage, distances, indices

    except Exception:
        traceback.print_exc()
        return 40.0, 100.0, 100.0,[],[]



