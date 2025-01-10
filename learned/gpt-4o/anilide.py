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
    
    # Define SMARTS pattern for anilide
    # Need an acyl group connected to an amine N, and that N is part of aniline structure
    # 'c1ccc(cc1)NC(=O)' - Phenyl ring connected to amide nitrogen
    anilide_smarts = Chem.MolFromSmarts("c1ccc(cc1)NC(=O)")
    
    # Check for anilide pattern in the given molecule
    if mol.HasSubstructMatch(anilide_smarts):
        return True, "Contains the distinct acylated aniline linkage of anilides"
    
    return False, "No distinct anilide motif detected"