"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:33569 anilide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline.

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
    
    # Look for benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"
    
    # Look for amide group attached to benzene ring
    amide_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group attached to benzene ring"
    
    # Look for acylated amide nitrogen
    acyl_pattern = Chem.MolFromSmarts("C(=O)NC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Amide nitrogen not acylated"
    
    # All conditions met, classify as anilide
    return True, "Contains aromatic amide group obtained by acylation of aniline"