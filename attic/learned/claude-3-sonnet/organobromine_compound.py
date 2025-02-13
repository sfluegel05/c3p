"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: organobromine compound
Definition: A compound containing at least one carbon-bromine bond.
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound must contain at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of both C and Br atoms
    has_carbon = False
    has_bromine = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            has_carbon = True
        elif atom.GetAtomicNum() == 35:  # Bromine
            has_bromine = True
        if has_carbon and has_bromine:
            break
    
    if not has_carbon:
        return False, "No carbon atoms found"
    if not has_bromine:
        return False, "No bromine atoms found"

    # Look for C-Br bond pattern
    cbr_pattern = Chem.MolFromSmarts("[#6]-[Br]")
    if not mol.HasSubstructMatch(cbr_pattern):
        return False, "No carbon-bromine bond found"

    # Count C-Br bonds
    cbr_matches = len(mol.GetSubstructMatches(cbr_pattern))
    
    return True, f"Found {cbr_matches} carbon-bromine bond{'s' if cbr_matches > 1 else ''}"