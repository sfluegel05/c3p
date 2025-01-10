"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins are characterized by a macrocyclic structure consisting of four 
    pyrrole rings connected mainly through methine bridges, forming a 16-membered macrocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # New SMARTS pattern for a macrocyclic structure typical of porphyrins.
    # The pattern varies to include potential porphyrin derivatives with different linkage styles.
    porphyrin_pattern = Chem.MolFromSmarts('n1ccccc1-[n+]-2C=CC=C-3C=CC=C[n+]3-2[nH]c4c[n+]ccc4')
    if porphyrin_pattern is None:
        return None, "Error in SMARTS pattern definition"

    # Attempt to find the macrocyclic structure
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "SMILES string contains a porphyrin-like macrocyclic structure"
    else:
        return False, "No macrocyclic structure typical of porphyrins found in SMILES string"