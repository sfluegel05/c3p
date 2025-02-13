"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is composed of a D-ribose sugar attached to a nucleobase via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Ribose sugar pattern (allowing for minor substitutions at OH groups)
    ribose_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O')
    # Nucleobase heterocycle pattern (purine or pyrimidine ring)
    nucleobase_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2')  # Purine base
    nucleobase_pattern_pyrimidine = Chem.MolFromSmarts('n1ccc(=O)[nH]c1=O')  # Pyrimidine base

    # Check for ribose sugar
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "No D-ribose sugar moiety found"

    # Check for nucleobase
    nucleobase_matches = mol.GetSubstructMatches(nucleobase_pattern)
    nucleobase_pyrimidine_matches = mol.GetSubstructMatches(nucleobase_pattern_pyrimidine)
    if not nucleobase_matches and not nucleobase_pyrimidine_matches:
        return False, "No nucleobase moiety found"

    # Check for glycosidic bond between sugar and base
    # Look for bonds between anomeric carbon of sugar and nitrogen of nucleobase
    anomeric_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4]
    base_nitrogens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetDegree() == 2]

    glycosidic_bond_found = False
    for ac in anomeric_carbons:
        for bn in base_nitrogens:
            bond = mol.GetBondBetweenAtoms(ac, bn)
            if bond:
                glycosidic_bond_found = True
                break
        if glycosidic_bond_found:
            break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond between sugar and base found"

    # Check for phosphate groups (to exclude nucleotides)
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)O')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if phosphate_matches:
        return False, "Phosphate group found; molecule is a nucleotide, not a nucleoside"

    return True, "Molecule is a ribonucleoside with D-ribose sugar and nucleobase connected via glycosidic bond"