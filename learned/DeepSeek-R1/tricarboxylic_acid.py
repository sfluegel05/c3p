"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:35711 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid has exactly three carboxy groups (-COOH or deprotonated form).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the carboxy group pattern (including both -COOH and -COO-)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    matches = mol.GetSubstructMatches(carboxy_pattern)
    
    if len(matches) == 3:
        # Verify that there are exactly three COOH/COO- groups and no more
        # Also check that there are no ester or amide groups which might indicate non-acidic COO
        ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
        amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        
        if len(ester_matches) + len(amide_matches) > 0:
            return False, "Found ester or amide groups, which may indicate non-carboxy COO"
        
        return True, "Contains exactly three carboxy groups"
    else:
        return False, f"Found {len(matches)} carboxy groups, expected 3"