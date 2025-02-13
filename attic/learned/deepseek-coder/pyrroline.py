"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:38943 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a dihydropyrrole ring, which is a five-membered ring with one nitrogen and one double bond.

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

    # Define the pyrroline pattern: a five-membered ring with one nitrogen and one double bond
    pyrroline_pattern = Chem.MolFromSmarts("[nH0]1[CH]=[CH][CH2]1")
    
    # Check if the molecule contains the pyrroline pattern
    if mol.HasSubstructMatch(pyrroline_pattern):
        return True, "Contains a dihydropyrrole ring (pyrroline)"
    else:
        return False, "No dihydropyrrole ring (pyrroline) found"