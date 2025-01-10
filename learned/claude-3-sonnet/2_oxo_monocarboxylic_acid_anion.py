"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:38478 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    Must have a carboxylate group with an oxo group at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-oxo carboxylate pattern
    # [O-]C(=O)C(=O)* pattern ensures the oxo group is exactly at position 2
    pattern = Chem.MolFromSmarts("[O-]C(=O)C(=O)[#6,#7,#8,#9,#16,#15,H]")
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No 2-oxo carboxylate group found"

    # For each match, verify it's a valid 2-oxo carboxylate
    for match in matches:
        carboxylate_carbon = match[1]  # C of C(=O)[O-]
        oxo_carbon = match[2]  # C of C(=O)R
        
        # Check that the oxo carbon is not part of another carboxylate
        # This prevents matching molecules where the pattern is part of other structures
        other_carboxylate = Chem.MolFromSmarts("[O-]C(=O)")
        other_matches = mol.GetSubstructMatches(other_carboxylate)
        is_valid = True
        
        for other_match in other_matches:
            if other_match[1] == oxo_carbon:  # If oxo carbon is part of another carboxylate
                is_valid = False
                break
        
        if is_valid:
            # Additional validation: ensure the oxo group is not part of an acid anhydride
            anhydride_pattern = Chem.MolFromSmarts("C(=O)OC(=O)")
            if not mol.HasSubstructMatch(anhydride_pattern):
                # Get the atom attached to the oxo group
                for neighbor in mol.GetAtomWithIdx(oxo_carbon).GetNeighbors():
                    if neighbor.GetIdx() not in [carboxylate_carbon]:
                        # If we found a valid 2-oxo carboxylate, return True
                        return True, "Contains carboxylate with oxo group at 2-position"

    return False, "No valid 2-oxo carboxylate group found"