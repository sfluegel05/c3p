"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

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

    # Remove hydrogens for simplicity
    mol = Chem.RemoveHs(mol)

    # Check for triterpenoid core (pentacyclic ring system)
    # Triterpenoids often have a characteristic pentacyclic structure
    pentacyclic_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCC5C4CCCC5')  # Simplified core
    if pentacyclic_pattern is None:
        return False, "Failed to parse triterpenoid core SMARTS pattern"

    if not mol.HasSubstructMatch(pentacyclic_pattern):
        return False, "No triterpenoid core detected"

    # Detect sugar moieties (glycosidic units)
    # Using a generic sugar pattern (pyranose ring)
    sugar_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C1O')
    if sugar_pattern is None:
        return False, "Failed to parse sugar SMARTS pattern"

    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moieties found"

    # Check for glycosidic linkage between core and sugar
    # Look for an oxygen atom connecting the core to the sugar
    core_atoms = mol.GetSubstructMatch(pentacyclic_pattern)
    sugar_atom_indices = [atom_idx for match in sugar_matches for atom_idx in match]
    core_atom_set = set(core_atoms)
    sugar_atom_set = set(sugar_atom_indices)

    found_linkage = False
    for bond in mol.GetBonds():
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Check for an oxygen linkage between core and sugar
        if bond.GetBondType() == Chem.BondType.SINGLE:
            if (atom1_idx in core_atom_set and atom2_idx in sugar_atom_set) or \
               (atom2_idx in core_atom_set and atom1_idx in sugar_atom_set):
                if atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8:
                    found_linkage = True
                    break

    if not found_linkage:
        return False, "No glycosidic linkage between triterpenoid core and sugar moiety found"

    # Optional: Check for additional hydroxyl or carboxyl groups common in saponins
    functional_groups = ['[OX2H]', 'C(=O)O']  # Hydroxyl and carboxyl groups
    has_functional_group = False
    for fg_smarts in functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if fg_pattern and mol.HasSubstructMatch(fg_pattern):
            has_functional_group = True
            break

    if not has_functional_group:
        return False, "No hydroxyl or carboxyl groups found; unusual for triterpenoid saponins"

    return True, "Molecule is a triterpenoid saponin with a triterpenoid core and glycosidically linked sugar moiety"