"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two, or three hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens for accurate valence calculation
    mol = Chem.AddHs(mol)
    
    # Get all nitrogen atoms
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogens:
        return False, "No nitrogen atoms present"
    
    # Preprocess exclusion patterns
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')  # Nitro groups
    amide_pattern = Chem.MolFromSmarts('C(=O)N')        # Amides
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    for atom in nitrogens:
        # Skip aromatic nitrogens (e.g., pyridine, pyrrole)
        if atom.GetIsAromatic():
            continue
        
        # Check if part of nitro group
        in_nitro = any(atom.GetIdx() in match for match in nitro_matches)
        if in_nitro:
            continue
        
        # Check if part of amide group
        in_amide = any(atom.GetIdx() in match for match in amide_matches)
        if in_amide:
            continue
        
        # Check neighbors: only carbon or hydrogen allowed
        carbon_neighbors = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in {1, 6}:  # Only H or C allowed
                break
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors.append(neighbor)
        else:  # Only executed if loop didn't break (all neighbors are H/C)
            # Verify number of carbon substituents (1-3)
            num_carbons = len(carbon_neighbors)
            if num_carbons not in {1, 2, 3}:
                continue
            
            # Check total valence (carbons + hydrogens) == 3
            num_h = atom.GetTotalNumHs()
            if (num_carbons + num_h) != 3:
                continue
            
            # Successfully found valid amine group
            amine_type = ['primary', 'secondary', 'tertiary'][num_carbons-1]
            return True, f"Contains {amine_type} amine group"
    
    return False, "No valid amine groups detected"