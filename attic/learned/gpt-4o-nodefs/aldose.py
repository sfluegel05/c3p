"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is characterized by an aldehyde group in an open chain form
    or can exist in cyclic forms such as pyranoses and furanoses.

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

    # Check for presence of aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")  # General aldehyde group
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Ensure presence of multiple hydroxyl groups
        hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) >= 3:
            return True, "Open-chain form with aldehyde and sufficient hydroxyls"
    
    # Check for six-membered cyclic (pyranose form) with oxygens and hydroxyls
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Cyclic pyranose form detected consistent with aldose"

    # Check for five-membered cyclic (furanose form)
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(CO)C1")
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Cyclic furanose form detected consistent with aldose"

    # Additional check: generic pattern for aldoses (open-chain or cyclic)
    generic_aldose_pattern = Chem.MolFromSmarts("C(O)(C(O)C(O)C=O)")
    if mol.HasSubstructMatch(generic_aldose_pattern):
        return True, "General aldose pattern matched"

    return False, "Structure does not fit typical aldose characteristics"