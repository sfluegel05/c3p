"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for 3-oxo group and 5alpha-steroid framework
    # 3-Oxo group (C=O) in the third position and key features of the steroid structure
    # The SMARTS pattern focuses on 3-oxo group presence and fused steroid rings with specific stereochemistry
    oxo_5alpha_pattern = Chem.MolFromSmarts("C1C2C3[C@H](CC[C@@H]3C(=O))C4CCC[C@@H]4C[C@@H]2CC1")
    
    if mol.HasSubstructMatch(oxo_5alpha_pattern):
        return True, "Matches the 3-oxo-5alpha-steroid pattern"
    else:
        return False, "Does not contain 3-oxo group in a 5alpha-steroid framework"

    return None, None