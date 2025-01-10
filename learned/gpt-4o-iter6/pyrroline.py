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

    # Refined pyrroline pattern: five-membered ring with specific nitrogen and double bond arrangements
    pyrroline_patterns = [
        "C1C=CN=C1",  # 1-Pyrroline with specific double bond
        "C1=CNCC1",   # 2-Pyrroline with specific bond positions
        "C1=CC=NC1",  # Allow aromatic pyrroline derivatives
        "C1=CCN=C1",  # Alt form capturing tautomeric shifts
        "C1=CN=CC1",  # Pyrrolone variations more definitively
        "[n;H1]1cccc1",  # Aromatic pyrroline recognition
        "n1cccc1",    # Non-hydrogen bonded variants for aromaticity check
    ]

    # Check for any pyrroline substructure match
    for pattern in pyrroline_patterns:
        pyrroline_mol = Chem.MolFromSmarts(pattern) # convert SMARTS pattern to molecule
        if pyrroline_mol and mol.HasSubstructMatch(pyrroline_mol):
            return True, f"Contains pyrroline-like structure: match for pattern {pattern}"

    return False, "No pyrroline-like structure detected"