import unittest
import numpy as np
from actseek.ActSeekLib_cy import  ransac_protein, Active_site
from actseek.ActSeek_cy import read_pdbs, printProtein, ActSeekMain, read_config 

class TestActSeek(unittest.TestCase):
    def test_ActSeek(self):
        aa_active_an=[292,448,478]
        aa_active=[0,1,2]
        seed_protein = "AF-Q9ZHI2-F1-model_v4.pdb"
        case_protein_name="P45370"
        case_protein_filename = "AF-P45370-F1-model_v4.pdb"
        protein_coords, protein_coords_cb, seed_coords,cavity_coords, cavity_coords_cb,  aa, aaCav, active, real_index, real_index_seed = read_pdbs(
            aa_active_an,
            seed_protein,
            case_protein_filename)

        index_used = np.array([int(x) for x in aa_active])
        active_used = np.array([active[x] for x in index_used])
        cavity_coords_used = np.array([cavity_coords[x] for x in index_used])
        cavity_coords_cb_used = np.array([cavity_coords_cb[x] for x in index_used])

        aa_des = {"GLY": "GLY", "ALA": "ALA", "PRO": "PRO", "ARG": "ARG", "HIS": "HIS", "LYS": "LYS", "ASP": "ASP",
                  "GLU": "GLU",
                  "SER": "SER", "THR": "MET", "ASN": "ASN", "GLN": "GLN", "CYS": "CYS", "VAL": "VAL", "ILE": "MET",
                  "LEU": "LEU",
                  "MET": "MET", "PHE": "PHE", "TYR": "TYR", "TRP": "TRP", "LIG": "DIF", "MSE": "HFBd", "HSD": "HIS"}
        active_site = Active_site(active_used, aaCav, cavity_coords_used, aa_des)
        active_site.get_AA_groups(False)

        pc = active_site.get_possible_correspondences(protein_coords, aa, real_index)
        threshold = 4
        pc2 = [[4, 1], [8, 1], [16, 1], [35, 1], [48, 1], [56, 1], [83, 1], [87, 1], [100, 1], [105, 1], [110, 1], [113, 1], [119, 1], [120, 1], [127, 1], [129, 0], [131, 1], [137, 2], [140, 1], [148, 0], [161, 2], [163, 1], [176, 1], [181, 1], [192, 1], [194, 1], [198, 1], [230, 1], [233, 1], [234, 1], [236, 1], [250, 1], [253, 1], [265, 1], [285, 1], [288, 1], [291, 0], [301, 1], [302, 2], [307, 1], [320, 1], [330, 2]]

        self.assertEqual(pc, pc2)

        valid_combinations = active_site.get_all_possible_combinations(pc, protein_coords, aa, real_index, threshold, aa_des)


        valid_conbinations2 = [([148, 0], [301, 1], [302, 2]), ([148, 0], [301, 1], [330, 2]), ([161, 2], [163, 1], [291, 0])]

        self.assertEqual(valid_combinations, valid_conbinations2)
        iterations=2000

        distances, distances_arround, t_transformed, solution, translation_vector, rotation = ransac_protein(
            cavity_coords, protein_coords, protein_coords_cb, seed_coords, valid_combinations,
            iterations, aa_des, aa, real_index, real_index_seed, aaCav, active, cavity_coords_used,cavity_coords_cb_used, 4,3, 0)

        self.assertEqual(distances, [np.float64(0.10003525945620756), np.float64(0.12675655227746951), np.float64(0.23646291635390213)])

        self.assertEqual(distances_arround, 0.9195135831832886 )     
        self.assertEqual(solution[0][0], 148)
        self.assertEqual(solution[0][1], 0)
        self.assertEqual(solution[1][0], 301)
        self.assertEqual(solution[1][1], 1)
        self.assertEqual(solution[2][0], 330)
        self.assertEqual(solution[2][1], 2)

        
        self.assertEqual(translation_vector[0], -4.844457852656266)
        self.assertEqual(translation_vector[1], -0.4984162025860206)
        self.assertEqual(translation_vector[2],  5.679361716090221)

    
        self.assertEqual(rotation[0][0], 0.6954747836157743)
        self.assertEqual(rotation[0][1], -0.3429827347761849)
        self.assertEqual(rotation[0][2], -0.6314092721840895)
        self.assertEqual(rotation[1][0], 0.7022705641664962)
        self.assertEqual(rotation[1][1], 0.5104380339521987)
        self.assertEqual(rotation[1][2], 0.4962550435011074)      
        self.assertEqual(rotation[2][0], 0.15208839554635092)
        self.assertEqual(rotation[2][1], -0.7885530147938461)
        self.assertEqual(rotation[2][2], 0.5958634598628055)


        self.assertEqual(t_transformed[0][0], 35.19056829964826)
        self.assertEqual(t_transformed[0][1], -17.331403019539664)
        self.assertEqual(t_transformed[0][2], -30.405356987032683)
        self.assertEqual(t_transformed[1][0], 33.04085405217916)
        self.assertEqual(t_transformed[1][1], -14.186041777374509)
        self.assertEqual(t_transformed[1][2], -30.88662358765822)
        self.assertEqual(t_transformed[2][0], 30.535122323687023)
        self.assertEqual(t_transformed[2][1], -14.628568829486099)
        self.assertEqual(t_transformed[2][2], -27.98975784447869)
        self.assertEqual(t_transformed[3][0], 27.41529875222943)
        self.assertEqual(t_transformed[3][1], -13.460922281360759)
        self.assertEqual(t_transformed[3][2], -29.94450059030705)
        self.assertEqual(t_transformed[148][0], 1.1990442794944398)
        self.assertEqual(t_transformed[148][1], 0.24718262365508192)
        self.assertEqual(t_transformed[148][2], 4.851651386124743)
        self.assertEqual(t_transformed[301][0], 10.378191029917668)
        self.assertEqual(t_transformed[301][1],-4.018548325546492)
        self.assertEqual(t_transformed[301][2], -0.1071972326115942)
        self.assertEqual(t_transformed[330][0], 5.83985406364274)
        self.assertEqual(t_transformed[330][1], -4.21989251978693)
        self.assertEqual(t_transformed[330][2], -1.3365764975664227)
        self.assertEqual(t_transformed[354][0], -7.092057499044784)
        self.assertEqual(t_transformed[354][1], -16.73018974792045)
        self.assertEqual(t_transformed[354][2], 19.488264657200983)

        
      

if __name__ == '__main__':
    unittest.main()
    

