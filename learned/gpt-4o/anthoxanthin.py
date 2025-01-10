"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with a flavonoid core and multiple hydroxyl, methoxy, or glycosyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an anthoxanthin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for flavonoid core (C6-C3-C6 structure)
    flavonoid_core_patterns = [
        Chem.MolFromSmarts('c1cc(-c2cccc(o2)c3=cc(=O)oc3)cc1'),  # Common flavone pattern
        Chem.MolFromSmarts('c1cc(-c2ccoc2)c(c(-c3ccccc3=O)c1)o'),  # Chalcone-like pattern
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_core_patterns):
        return False, "No recognizable flavonoid core structure found"
    
    # Search for the presence of at least two hydroxyl groups (common in flavonoids)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyls < 2:
        return False, f"Insufficient number of hydroxy groups, found {num_hydroxyls}"

    # Search for methoxy groups (common substitution in flavonoids)
    methoxy_pattern = Chem.MolFromSmarts('CO')
    num_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))
    
    # Search for possible glycoside linkages (attachment of sugar residues)
    glycoside_pattern = Chem.MolFromSmarts('[C@H]1(O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O)')
    num_glycosides = len(mol.GetSubstructMatches(glycoside_pattern))
    
    # Consider the anthoxanthins criteria
    if num_glycosides > 0 or num_methoxy > 0:
        return True, "Contains recognizable flavonoid core with hydroxy, methoxy or glycoside groups, indicating potential anthoxanthin structure"

    # Evaluate based on the variety of functional groups
    if num_hydroxyls >= 2 and (num_methoxy > 0 or num_glycosides > 0):
        return True, "Contains flavonoid core with diverse functional groups typical of anthoxanthins"
    
    return False, "Lacks sufficient structure or functional group diversity for classification as anthoxanthin"