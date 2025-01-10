"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is an organic heteromonocyclic compound with a structure based on a dihydropyrrole 
    (a five-membered ring containing one nitrogen and at least one double bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pyrroline pattern: five-membered ring with one nitrogen and double bond(s)
    pyrroline_pattern = Chem.MolFromSmarts("C1CC=NC1")  # 1-pyrroline as a template

    # Check for a pyrroline substructure match
    if mol.HasSubstructMatch(pyrroline_pattern):
        return True, "Contains pyrroline structure (five-membered ring with nitrogen and double bond)"
    else:
        return False, "No pyrroline structure detected"