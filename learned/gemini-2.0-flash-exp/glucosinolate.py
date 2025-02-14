"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # 1. Check for thioglucoside core (pyranose linked to sulfur)
    thioglucoside_pattern = Chem.MolFromSmarts("[C]1[O][C][C][C][C]1[S]")
    thioglucoside_matches = mol.GetSubstructMatches(thioglucoside_pattern)
    if not thioglucoside_matches:
        return False, "Thioglucoside core not found"

    # 2. Check for the central C atom bonded to S, N, and a side-chain
    central_carbon_pattern = Chem.MolFromSmarts("[S][C](=[N])")
    central_carbon_matches = mol.GetSubstructMatches(central_carbon_pattern)
    if not central_carbon_matches:
        return False, "Central carbon with S and N not found"

    # 3. Check for sulfonated oxime group (anionic and non-anionic)
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[C]=[N][O]S(=O)(=O)[O-,O]")
    sulfonated_oxime_matches = mol.GetSubstructMatches(sulfonated_oxime_pattern)
    if not sulfonated_oxime_matches:
      return False, "Sulfonated oxime group not found"

    # 4. Check connectivity between thioglucoside S and central C
    
    central_c_atoms = [match[1] for match in central_carbon_matches] # indices of the central carbons
    sulfur_atoms_central_group = [mol.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx() for idx in central_c_atoms] #indices of the sulfur atoms directly connected to the central carbon
    thioglucoside_s_atoms = [match[-1] for match in thioglucoside_matches] # indices of the S atoms of the thioglucoside

    if not any(s_center in thioglucoside_s_atoms for s_center in sulfur_atoms_central_group):
        return False, "Thioglucoside sulfur not connected to the central carbon"
    
    #5. Check for C=N bond stereochemistry

    for match in sulfonated_oxime_matches:
        c_idx = match[0]
        n_idx = match[1]
        o_idx = match[2]
    
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)

        c_neighbors = [neighbor.GetIdx() for neighbor in c_atom.GetNeighbors()]
        n_neighbors = [neighbor.GetIdx() for neighbor in n_atom.GetNeighbors()]
            
        #Find neighbors of C that are not N or H. If there are no neighbors return False because the molecule doesn't have an anti configuration.
        c_other_neighbor_index = None
        for c_neigh_index in c_neighbors:
            if c_neigh_index != n_idx and mol.GetAtomWithIdx(c_neigh_index).GetAtomicNum() != 1:
                c_other_neighbor_index = c_neigh_index
                break

        if c_other_neighbor_index is None:
           return False, "No anti configuration found on C=N bond"

        n_other_neighbor_index = None
        for n_neigh_index in n_neighbors:
            if n_neigh_index != o_idx and mol.GetAtomWithIdx(n_neigh_index).GetAtomicNum() != 1:
                n_other_neighbor_index = n_neigh_index
                break

        if n_other_neighbor_index is None:
           return False, "No anti configuration found on C=N bond"

        # check anti configuration using bond direction
        c_neigh_pos = mol.GetConformer().GetAtomPosition(c_other_neighbor_index)
        n_neigh_pos = mol.GetConformer().GetAtomPosition(n_other_neighbor_index)
        c_pos = mol.GetConformer().GetAtomPosition(c_idx)
        n_pos = mol.GetConformer().GetAtomPosition(n_idx)

        vec1 = c_neigh_pos - c_pos
        vec2 = n_neigh_pos - n_pos
        cross_product = vec1[0] * vec2[1] - vec1[1] * vec2[0]
        
        if cross_product > 0:
            return False, "Anti configuration not found on C=N bond"
        
    return True, "Molecule matches all required patterns for glucosinolate"