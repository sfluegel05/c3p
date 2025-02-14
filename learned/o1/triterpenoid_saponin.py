"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside where the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common triterpenoid skeletons
    triterpene_patterns = [
        # Oleanane skeleton
        Chem.MolFromSmarts('C1CCC2(C)C3CCC4(C)C(CCC5C(C)(C)CCC(=O)C5(C)C)C4C3CCC12C'),
        # Ursane skeleton
        Chem.MolFromSmarts('C1CCC2(C)C3CCC4(C)C(CCC5C(C)(C)CCC(=O)C5(C)C)C4C3CCC12C'),
        # Lupane skeleton
        Chem.MolFromSmarts('C1CCC2(C)C3C4CCC5(C)C(CCC6C(C)(C)CCC(=O)C6(C)C)C5C4CCC3C12'),
        # Other triterpenoid skeletons can be added here
    ]

    # Check for triterpenoid core
    triterpenoid_core_found = False
    for pattern in triterpene_patterns:
        if mol.HasSubstructMatch(pattern):
            triterpenoid_core_found = True
            break

    if not triterpenoid_core_found:
        return False, "No triterpenoid core found"

    # Define SMARTS pattern for sugar moieties (glycosidic units)
    sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O')  # Generic pyranose ring

    # Find all sugar moieties
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties found"

    # Check for glycosidic bond between triterpenoid core and sugar
    # Find all oxygen atoms connected to the core and sugars
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]

    glycosidic_bond_found = False
    for oxy_atom in oxygen_atoms:
        neighbors = oxy_atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        atom1, atom2 = neighbors
        atom1_in_core = any(mol.HasSubstructMatch(pattern, atoms=[atom1.GetIdx()]) for pattern in triterpene_patterns)
        atom2_in_sugar = mol.HasSubstructMatch(sugar_pattern, atoms=[atom2.GetIdx()])
        if atom1_in_core and atom2_in_sugar or atom2_in_core and atom1_in_sugar:
            glycosidic_bond_found = True
            break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond found between triterpenoid core and sugar moiety"

    return True, "Contains triterpenoid core with glycosidically linked sugar moiety"