"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is characterized by a five-membered ring containing a nitrogen atom and at least one C=C double bond.

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

    # Define a general pyrroline SMARTS pattern: a non-aromatic five-membered ring with a nitrogen and at least one double bond
    pyrroline_patterns = [
        Chem.MolFromSmarts("C1=NC=CC1"),  # 2-pyrroline
        Chem.MolFromSmarts("C1=CC=NC1"),  # 1-pyrroline
        Chem.MolFromSmarts("C1=CCC=N1")   # 3-pyrroline
    ]
    
    # Check if the molecule contains any pyrroline substructure
    for pattern in pyrroline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a pyrroline structure: five-membered ring with a nitrogen and a double bond"
    
    return False, "Lacks pyrroline structure: specific five-membered ring criteria not satisfied"