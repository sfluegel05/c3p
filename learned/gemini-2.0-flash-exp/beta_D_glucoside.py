"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a D-glucopyranose ring
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([CH2]O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)O1")
    matches = mol.GetSubstructMatches(glucose_pattern)
    if not matches:
       return False, "No D-glucose ring found."

    for match in matches:
        # Get anomeric carbon (C1) index
        anomeric_carbon_index = match[0]
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_index)

        #check that anomeric carbon is a CH group, otherwise return False
        if anomeric_carbon.GetTotalNumHs() != 1:
            continue
        
        # find the two oxygen neighbours of the anomeric carbon.
        oxygen_neighbors = []
        for neighbor in anomeric_carbon.GetNeighbors():
             if neighbor.GetAtomicNum() == 8:
                oxygen_neighbors.append(neighbor)
        
        if len(oxygen_neighbors) != 2:
            continue
        
        # find the ring oxygen, which is the one attached to the ring, and the other (glycosidic) oxygen
        ring_oxygen = None
        glycosidic_oxygen = None
        for oxygen in oxygen_neighbors:
            if oxygen.GetIdx() in match:
                ring_oxygen = oxygen
            else:
                glycosidic_oxygen = oxygen
        
        if ring_oxygen is None or glycosidic_oxygen is None:
          continue
            
        # get indices of atoms attached to the anomeric carbon.
        anomeric_neighbors_indices = [n.GetIdx() for n in anomeric_carbon.GetNeighbors()]
        
        # find the index of the H at the anomeric carbon
        for i in anomeric_neighbors_indices:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() == 1: # found hydrogen!
                h_index = i
                break
        else:
            continue #should not happen
        
        # find the C5 carbon, which is the last carbon in the ring (not CH2OH)
        c5_index = match[4]
        c5_atom = mol.GetAtomWithIdx(c5_index)

        # find the CH2OH carbon at C6, which is bonded to C5
        c6_index = None
        for neighbor in c5_atom.GetNeighbors():
          if neighbor.GetIdx() not in match and neighbor.GetAtomicNum() == 6:
             c6_index = neighbor.GetIdx()
             break
        else:
            continue

        # Find one of the H at C6
        c6_atom = mol.GetAtomWithIdx(c6_index)
        for neighbor in c6_atom.GetNeighbors():
          if neighbor.GetAtomicNum() == 1:
              h6_index = neighbor.GetIdx()
              break
        else:
             continue

        # Get the position of the atoms to check for stereochemistry
        anomeric_carbon_coords = mol.GetConformer().GetAtomPosition(anomeric_carbon_index)
        h_coords = mol.GetConformer().GetAtomPosition(h_index)
        c6_coords = mol.GetConformer().GetAtomPosition(c6_index)
        h6_coords = mol.GetConformer().GetAtomPosition(h6_index)

        # Calculate vectors
        vec_h = h_coords - anomeric_carbon_coords
        vec_ch2oh = h6_coords - c6_coords
        
        # project vec_h onto plane defined by the C1-C2-C5-C6-h atoms.
        plane_normal = (mol.GetConformer().GetAtomPosition(match[1]) - anomeric_carbon_coords).CrossProduct(mol.GetConformer().GetAtomPosition(c5_index) - anomeric_carbon_coords) # C2-C1 cross C5-C1
        proj_h = vec_h - (vec_h.DotProduct(plane_normal) / plane_normal.LengthSq()) * plane_normal

        # project vec_ch2oh onto the same plane
        proj_ch2oh = vec_ch2oh - (vec_ch2oh.DotProduct(plane_normal) / plane_normal.LengthSq()) * plane_normal

        # Check the dot product to verify if they are on the same side.
        dot_prod = proj_h.DotProduct(proj_ch2oh)
        if dot_prod > 0: #H and CH2OH are on the same side, thus beta configuration
            return True, "beta-D-glucoside detected (glycosidic bond with correct beta configuration)"

    return False, "Not a beta-D-glucoside (no glycosidic bond with correct beta configuration)"