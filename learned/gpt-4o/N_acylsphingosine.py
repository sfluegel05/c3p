"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is composed of sphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update sphingosine backbone pattern:
    # Looking for a long chain structure with key functional groups
    sphingosine_pattern = Chem.MolFromSmarts("NC(CO)C(O)CCCCCCCC=C")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Robust check for amide bond attached to nitrogen
    amide_pattern = Chem.MolFromSmarts("N[C](=O)[C]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found connecting nitrogen to fatty acyl group"

    return True, "Contains sphingosine backbone with an amide bond connecting to a fatty acyl group"