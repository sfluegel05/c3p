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

    # Modified aldehyde pattern to account for potential hemiacetal forms
    aldehyde_or_hemiacetal_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6] | [O][CX4H2][CX3](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_or_hemiacetal_pattern):
        return False, "No aldehyde or cyclic hemiacetal forms detected"

    # Check for multiple hydroxyl groups (polyhydroxy structure)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient number of hydroxyl groups for polyhydroxy structure"

    # Checking for cyclic forms - furanose or pyranose
    furanose_or_pyranose_pattern = Chem.MolFromSmarts("[O][C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1 | [O][C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1")
    cyclic_match = mol.HasSubstructMatch(furanose_or_pyranose_pattern)
    if not cyclic_match:
        return False, "No matching 5 or 6-membered cyclic forms typical of aldoses"

    return True, "Contains structural features consistent with aldose, including polyhydroxy and potential cyclic (hemiacetal) forms"