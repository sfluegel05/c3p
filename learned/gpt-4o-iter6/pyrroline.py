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
        return None, "Invalid or non-parsable SMILES string"

    # Improved pyrroline pattern: focusing on five-membered rings with nitrogen and double bond
    pyrroline_patterns = [
        "C1=CC=CN1",  # Pyrrole-like structures with flexibility in double bonds
        "C1C=CCN1",   # Ensuring the nitrogen and one double bond
        "C1=CCNC1",   # Accounting for alternative double bond positions
        "N1=CC=CC1",  # Including different tautomeric forms
        "C1=NC=CC1",  # Another ring formation pattern
    ]

    # Check all substructure patterns
    for pattern in pyrroline_patterns:
        pyrroline_mol = Chem.MolFromSmarts(pattern)
        if pyrroline_mol and mol.HasSubstructMatch(pyrroline_mol):
            return True, f"Contains pyrroline-like structure: matched pattern {pattern}"

    return False, "No pyrroline structure detected"