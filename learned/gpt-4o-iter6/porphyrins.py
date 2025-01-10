"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is characterized by a macrocyclic structure of four pyrrole nuclei 
    connected through the alpha-positions, often but not exclusively, by four methine bridges.

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
    
    # SMARTS pattern for a basic porphyrin-like macrocyclic structure with pyrroles
    # Allows connectivity variances commonly seen in porphyrin derivatives
    porphyrin_pattern = Chem.MolFromSmarts('n1ccc(c2ccc[nH]2)c2[nH]ccc(c3ccc[nH]3)c2c1')
    if porphyrin_pattern is None:
        return None, "Error in SMARTS pattern definition"

    # Attempt to find the macrocyclic structure
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "SMILES string contains a porphyrin-like macrocyclic structure with interconnected pyrroles"
    else:
        return False, "No macrocyclic structure typical of porphyrins found in SMILES string"