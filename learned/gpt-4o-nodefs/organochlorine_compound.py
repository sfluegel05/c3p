"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains carbon and at least one chlorine atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of chlorine atoms
    has_chlorine = any(atom.GetAtomicNum() == 17 for atom in mol.GetAtoms())
    
    if has_chlorine:
        return True, "Contains chlorine atoms, classified as an organochlorine compound"
    else:
        return False, "No chlorine atoms found, not an organochlorine compound"