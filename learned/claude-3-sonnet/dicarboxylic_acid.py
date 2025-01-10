"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains exactly two carboxyl (-COOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group pattern (-C(=O)OH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Count number of carboxyl groups
    num_carboxyl = len(matches)
    
    if num_carboxyl == 0:
        return False, "No carboxyl groups found"
    elif num_carboxyl == 1:
        return False, "Only one carboxyl group found"
    elif num_carboxyl > 2:
        return False, f"Found {num_carboxyl} carboxyl groups, expected exactly 2"
    
    # Additional check to ensure the carboxyl groups are not part of the same carboxy anhydride group
    # Carboxy anhydrides have pattern -C(=O)OC(=O)-
    anhydride_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX3](=[OX1])")
    if mol.HasSubstructMatch(anhydride_pattern):
        return False, "Contains carboxy anhydride group instead of separate carboxyl groups"
    
    # Additional check to ensure the carboxyl groups are not part of ester linkages
    # Look for non-hydrogen atoms connected to the oxygen of the carboxyl group
    for match in matches:
        c_idx = match[0]  # Index of the carbon atom
        o_idx = match[2]  # Index of the hydroxyl oxygen atom
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Check if oxygen is connected to anything other than H and the carboxyl carbon
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetIdx() != c_idx and neighbor.GetAtomicNum() != 1:
                return False, "One or more carboxyl groups are part of ester linkages"

    return True, "Contains exactly two free carboxyl groups"