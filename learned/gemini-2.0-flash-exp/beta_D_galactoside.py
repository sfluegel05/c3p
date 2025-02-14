"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside contains a D-galactose ring with beta-configuration at its anomeric center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the D-galactose pyranose ring.
    # Specifically look for a 6 membered ring with the correct stereochemistry and OH groups in the right configuration.
    # Note, that the bond order between the ring O and anomeric C is unspecified to accommodate different glycosidic bonds (O, C).
    galactose_ring_pattern = Chem.MolFromSmarts("[#6][#8X2]1[#6H][#6@H](O)[#6@@H](O)[#6@@H]([#6H]2[#8H][#6H])1") # this matches D-galactose but has an unspecified anomeric carbon

    if not mol.HasSubstructMatch(galactose_ring_pattern):
       return False, "No D-galactose ring found"

    matches = mol.GetSubstructMatches(galactose_ring_pattern)

    for match in matches:
        # Extract the relevant atoms indices
        anomeric_carbon = match[0] #The anomeric carbon
        ring_oxygen = match[1] #Ring oxygen
        carbon2 = match[2]
        carbon5 = match[5]
        
        # Check if the anomeric carbon is indeed connected to the ring oxygen
        bond = mol.GetBondBetweenAtoms(anomeric_carbon, ring_oxygen)
        if bond is None:
          continue #This match is not valid, the carbon and oxygen might be somewhere else, check next match

        # Check if the stereochemistry at the anomeric carbon is correct (beta)
        # Beta is cis to CH2OH on C5 and trans to -OH at C2
        c2_neighbor_atoms = [x.GetIdx() for x in mol.GetAtom(carbon2).GetNeighbors() if x.GetIdx() not in match]
        if len(c2_neighbor_atoms) == 0:
           continue
        c2_neighbor_atom = c2_neighbor_atoms[0]

        c5_neighbor_atoms = [x.GetIdx() for x in mol.GetAtom(carbon5).GetNeighbors() if x.GetIdx() not in match and x.GetAtomicNum() == 6]
        if len(c5_neighbor_atoms) == 0:
           continue
        c5_neighbor_atom = c5_neighbor_atoms[0]

        #Check stereochemistry at the anomeric carbon
        #This can be done by computing the dihedral angle
        conf = mol.GetConformer()
        anomeric_carbon_pos = conf.GetAtomPosition(anomeric_carbon)
        oxygen_pos = conf.GetAtomPosition(ring_oxygen)
        c2_pos = conf.GetAtomPosition(carbon2)
        c2_neighbor_pos = conf.GetAtomPosition(c2_neighbor_atom)
        c5_pos = conf.GetAtomPosition(carbon5)
        c5_neighbor_pos = conf.GetAtomPosition(c5_neighbor_atom)
            
        #Get direction vector between atoms
        v1 = (oxygen_pos.x - anomeric_carbon_pos.x, oxygen_pos.y - anomeric_carbon_pos.y, oxygen_pos.z - anomeric_carbon_pos.z)
        v2 = (c2_pos.x - anomeric_carbon_pos.x, c2_pos.y - anomeric_carbon_pos.y, c2_pos.z - anomeric_carbon_pos.z)
        v3 = (c2_neighbor_pos.x - c2_pos.x, c2_neighbor_pos.y - c2_pos.y, c2_neighbor_pos.z - c2_pos.z)
        v4 = (c5_pos.x - anomeric_carbon_pos.x, c5_pos.y - anomeric_carbon_pos.y, c5_pos.z - anomeric_carbon_pos.z)
        v5 = (c5_neighbor_pos.x - c5_pos.x, c5_neighbor_pos.y - c5_pos.y, c5_neighbor_pos.z - c5_pos.z)

        import numpy as np

        def angle(v1,v2):
            v1_norm = np.linalg.norm(v1)
            v2_norm = np.linalg.norm(v2)
            dot = np.dot(v1, v2)
            if v1_norm == 0 or v2_norm == 0:
                return None
            cos_angle = dot / (v1_norm * v2_norm)
            return np.degrees(np.arccos(cos_angle))

        ang_c1_c2 = angle(v1,v2)
        ang_c1_c5 = angle(v1,v4)
        ang_c2 = angle(v2,v3)
        ang_c5 = angle(v4,v5)
    
        if ang_c1_c2 is None or ang_c1_c5 is None or ang_c2 is None or ang_c5 is None:
            continue
        
        #Check for the anti relationship of the H at C1 vs the hydroxyl at C2
        if ang_c2 > 110 and ang_c1_c2 < 70:
            if ang_c5 > 110 and ang_c1_c5 < 70:
                return True, "Contains a beta-D-galactoside moiety"

    return False, "Not a beta-D-galactoside, based on anomeric carbon stereochemistry"