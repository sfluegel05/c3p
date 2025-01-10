"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with a flavonoid core and varying functional groups including hydroxyls, methoxyls, and glycosides.

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

    # Expanded set of flavonoid core patterns
    flavonoid_core_patterns = [
        Chem.MolFromSmarts('c1cc(-c2cccc(o2)c3=cc(=O)oc3)cc1'),  # Common flavone pattern
        Chem.MolFromSmarts('c1cc(-c2ccoc2)c(c(-c3ccccc3=O)c1)o'),  # Chalcone-like pattern
        Chem.MolFromSmarts('c1c2c(cc3ccccc3)ooc2cc(=O)o1'),  # Isoflavone-like pattern
    ]
    
    # Check for flavonoid core
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_core_patterns):
        return False, "No recognizable flavonoid core structure found"
    
    # Expanded search for hydroxyl groups (OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyls < 2:
        return False, f"Insufficient number of hydroxyl groups, found {num_hydroxyls}"

    # Search for methoxy groups (CO)
    methoxy_pattern = Chem.MolFromSmarts('CO')
    num_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))

    # Search for glycosidic patterns with flexibility in patterns
    glycoside_pattern = Chem.MolFromSmarts('[C@H]1(O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O)')
    num_glycosides = len(mol.GetSubstructMatches(glycoside_pattern))

    # Analyze carbon and oxygen for ring formations indicative of glycosides
    num_carbon_oxygen_rings = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in {6, 8} and len(atom.GetNeighbors()) > 2)

    # Evaluate potential anthoxanthin based on group diversity and ring structures
    if (num_glycosides > 0 or num_methoxy >= 1) and num_hydroxyls >= 2:
        return True, "Recognizable flavonoid core with hydroxyl, methoxy or glycoside groups typical of anthoxanthins"

    return False, "Lacks sufficient structure or functional group diversity for classification as anthoxanthin"