"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:32952 primary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule contains a primary amine group based on its SMILES string.
    A primary amine has an NH2 group connected to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to properly count them on N atoms
    mol = Chem.AddHs(mol)

    # Primary amine pattern: N with 2 H atoms connected to C
    # Exclude cases where N is part of amide (N-C=O)
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][CX4]")
    amide_pattern = Chem.MolFromSmarts("[NX3H2][CX3](=[OX1])")
    
    # Find matches
    amine_matches = mol.GetSubstructMatches(primary_amine_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Convert matches to sets of nitrogen atoms for comparison
    amine_nitrogens = {match[0] for match in amine_matches}
    amide_nitrogens = {match[0] for match in amide_matches}
    
    # Remove amide nitrogens from amine matches
    true_primary_amines = amine_nitrogens - amide_nitrogens
    
    if not true_primary_amines:
        return False, "No primary amine groups found"
        
    # Additional checks for each potential primary amine nitrogen
    for n_idx in true_primary_amines:
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check number of hydrogens
        if n_atom.GetTotalNumHs() != 2:
            continue
            
        # Check number of heavy atom neighbors (should be exactly 1)
        if len(n_atom.GetNeighbors()) != 1:
            continue
            
        # Check if the neighbor is carbon
        neighbor = n_atom.GetNeighbors()[0]
        if neighbor.GetAtomicNum() != 6:  # 6 is atomic number for carbon
            continue
            
        # If we get here, we've found a valid primary amine
        return True, f"Contains primary amine group(s)"
        
    return False, "No valid primary amine groups found"