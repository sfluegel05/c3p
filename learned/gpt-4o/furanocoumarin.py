"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furochromene that consists of a furan ring fused with a coumarin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccco1")
    # Define SMARTS pattern for the coumarin structure
    coumarin_pattern = Chem.MolFromSmarts("c1cc2ccoc2cc1=O")
    
    # Check for furan and coumarin rings
    has_furan = mol.HasSubstructMatch(furan_pattern)
    has_coumarin = mol.HasSubstructMatch(coumarin_pattern)
    
    if not (has_furan and has_coumarin):
        return False, "The structure lacks necessary furan and/or coumarin moieties"
    
    # Check for fusion by looking for overlapping or contiguous furan and coumarin components
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    
    # Simplified approach: check if any atom from furan matches overlaps with coumarin matches
    for furan_atoms in furan_matches:
        for coumarin_atoms in coumarin_matches:
            fusion_atoms = set(furan_atoms).intersection(set(coumarin_atoms))
            if len(fusion_atoms) >= 1:  # At least one atom must be in common for fusion
                return True, "Contains fused furan and coumarin rings"
    
    return False, "Furan and coumarin rings are not fused"