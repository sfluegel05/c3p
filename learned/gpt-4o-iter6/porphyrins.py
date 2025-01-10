"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins are natural pigments characterized by a macrocyclic structure
    consisting of four pyrrole rings connected through methine bridges 
    (-CH=) forming a large conjugated 16-membered macrocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a new, broader SMARTS pattern for porphyrins
    # n: pyrrole nitrogen, c: pyrrole carbon, C: methine (-CH=)
    porphyrin_pattern = Chem.MolFromSmarts('n1ccc2[nH]c(c1)C=Cc3c[nH]c(c3)C=Cc2C=Cc4c[nH]cc4')

    if porphyrin_pattern is None:
        return None, "Error in SMARTS pattern definition"

    # Attempt to find the macrocyclic structure
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "SMILES string contains a porphyrin-like macrocyclic structure"
    else:
        return False, "No macrocyclic structure typical of porphyrins found in SMILES string"