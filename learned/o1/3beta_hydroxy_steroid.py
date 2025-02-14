"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36804 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group at the 3-position in the beta orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure stereochemistry is assigned
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Steroid nucleus SMARTS pattern (cyclopentanoperhydrophenanthrene core)
    steroid_pattern = Chem.MolFromSmarts('C1CC2CCC3C(C2C1)CC4CCC(C3)C4')  # Simplified steroid backbone
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid core not found"

    # Define the 3beta-hydroxy group pattern
    # Looking for a beta-oriented hydroxy group at position 3
    beta_oh_pattern = Chem.MolFromSmarts('[C@@H](O)[C@H]')  # beta-hydroxy at chiral carbon
    if beta_oh_pattern is None:
        return False, "Invalid 3beta-hydroxy SMARTS pattern"

    # Find all matches for the pattern
    matches = mol.GetSubstructMatches(beta_oh_pattern, useChirality=True)
    if not matches:
        return False, "No 3beta-hydroxy group found"

    # Verify that the hydroxy group is at the 3-position of the steroid nucleus
    for match in matches:
        hydroxyl_carbon_idx = match[0]
        # Check if the carbon is part of the steroid core
        if mol.GetAtomWithIdx(hydroxyl_carbon_idx).IsInRing():
            # Assuming that the numbering starts from ring A
            return True, "Contains steroid backbone with 3beta-hydroxy group"

    return False, "3beta-hydroxy group not at correct position in steroid core"