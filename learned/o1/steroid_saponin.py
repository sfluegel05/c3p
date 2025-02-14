"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a compound consisting of a hydroxysteroid backbone with sugar moieties attached via glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define steroid nucleus pattern (cyclopentanoperhydrophenanthrene nucleus)
    steroid_smarts = 'C1CCC2C1(CC[C@@H]3CC[C@H]4[C@H]3CC[C@@H]4C2)'

    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl groups attached to the steroid backbone
    hydroxyl_smarts = '[#6]-[CH]-[OH]'
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    if hydroxyl_pattern is None:
        return False, "Invalid hydroxyl SMARTS pattern"
    
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found on the steroid backbone"

    # Define sugar (glycoside) pattern (pyranose ring)
    sugar_smarts = 'C1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)'
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if sugar_pattern is None:
        return False, "Invalid sugar SMARTS pattern"

    # Check for sugar moieties attached via glycosidic bonds
    glycosidic_bond_smarts = '[C;!R][O][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    glycosidic_bond_pattern = Chem.MolFromSmarts(glycosidic_bond_smarts)
    if glycosidic_bond_pattern is None:
        return False, "Invalid glycosidic bond SMARTS pattern"

    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No sugar moieties attached via glycosidic bonds found"

    # Ensure that the sugar is attached to the steroid backbone
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    connected = False
    
    for s_match in steroid_matches:
        s_atoms = set(s_match)
        for sugar_match in sugar_matches:
            sugar_atoms = set(sugar_match)
            # Check if any bond connects the steroid and sugar atoms
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if (begin_idx in s_atoms and end_idx in sugar_atoms) or (begin_idx in sugar_atoms and end_idx in s_atoms):
                    connected = True
                    break
            if connected:
                break
        if connected:
            break
    if not connected:
        return False, "Steroid backbone not connected to sugar moiety via glycosidic bond"

    return True, "Contains a hydroxysteroid backbone with sugar moieties attached via glycosidic bonds"