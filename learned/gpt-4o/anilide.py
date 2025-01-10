"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as an aromatic amide obtained by acylation of aniline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for an aromatic amide
    # This is a carbonyl connected to nitrogen and the nitrogen should be part of an aromatic system
    anilide_smarts = Chem.MolFromSmarts("aN-C(=O)-")
    
    # Find the anilide pattern
    if mol.HasSubstructMatch(anilide_smarts):
        return True, "Contains an aromatic amide linkage characteristic of anilides"
    
    return False, "No anilide linkage detected"