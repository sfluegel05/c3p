"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # SMARTS pattern for any nitrogen atom to capture polyamines
    nitrogen_pattern = Chem.MolFromSmarts("[NX3,NX4+]")
    
    # Identifying all nitrogen matches
    nitrogen_matches = mol.GetSubstructMatches(nitrogen_pattern)
    
    # Count number of nitrogen atoms
    num_nitrogen_atoms = len(nitrogen_matches)
    
    if num_nitrogen_atoms >= 2:
        return True, f"Molecule contains {num_nitrogen_atoms} nitrogen atoms, indicating potential polyamine structure."
    else:
        return False, f"Molecule contains {num_nitrogen_atoms} nitrogen atom(s), fewer than required for a polyamine."