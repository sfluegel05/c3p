"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a polyhydroxy aldehyde or its intramolecular hemiacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify simple aldehyde group or hemiacetal cyclic structure
    # Linear aldehyde pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    # Cyclic hemiacetal pattern - generalized to 5 or 6-membered rings
    cyclic_pattern = Chem.MolFromSmarts("[O]1[C@H](O)[C@H](O)[C@H](O)[C@H]1 | [O]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1")
    
    # Check for either aldehyde or hemiacetal cyclic structure
    if not (mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(cyclic_pattern)): 
        return False, "No aldehyde group or hemiacetal cyclic form detected"

    # Check for multiple hydroxyl groups (polyhydroxy structure)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient number of hydroxyl groups for polyhydroxy structure"

    return True, "Contains structural features consistent with aldose, either polyhydroxy linear aldehyde or cyclic (hemiacetal) form"