"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is composed of a D-ribose sugar (possibly modified) attached to a nucleobase via a glycosidic bond.
    The sugar component must be D-ribose, which is a five-membered ring with four carbons and one oxygen,
    and specific stereochemistry corresponding to D-ribose.

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
    # Flexible D-ribose pattern allowing substitutions at hydroxyl groups
    ribose_pattern = Chem.MolFromSmarts("""
    [C@H]1([#6,#7,#8,#9])[C@@H]([#6,#7,#8,#9])[C@H]([#6,#7,#8,#9])[C@@H](CO)[O]1
    """)
    # Purine base pattern (including modifications)
    purine_pattern = Chem.MolFromSmarts('n1c[nH,c]c2c1n[c,n][c,n]c2')
    # Pyrimidine base pattern (including modifications)
    pyrimidine_pattern = Chem.MolFromSmarts('c1c([nH,c])[nH,c]c([nH,c])c1=O')

    # Check for ribose sugar ring with correct stereochemistry
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No D-ribose sugar moiety found or incorrect stereochemistry"

    # Check for nucleobase (purine or pyrimidine)
    if not mol.HasSubstructMatch(purine_pattern) and not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No nucleobase moiety found"

    # Exclude molecules with phosphate groups (nucleotides)
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)O')
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group found; molecule is a nucleotide, not a nucleoside"

    # Check for glycosidic bond between C1' of ribose and N of nucleobase
    # Find the anomeric carbon (C1') in ribose
    ribose_match = mol.GetSubstructMatch(ribose_pattern)
    if not ribose_match:
        return False, "No ribose sugar moiety found"

    anomeric_carbon = None
    for idx in ribose_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            # Anomeric carbon is connected to the ring oxygen and a substituent (nucleobase)
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
            neighbor_ats = [mol.GetAtomWithIdx(nbr_idx).GetAtomicNum() for nbr_idx in neighbors]
            if 8 in neighbor_ats and len(neighbors) > 2:
                anomeric_carbon = atom
                break
    if anomeric_carbon is None:
        return False, "Anomeric carbon not found in ribose ring"

    # Check if anomeric carbon is connected to a nitrogen in the nucleobase
    glycosidic_bond_found = False
    for bond in anomeric_carbon.GetBonds():
        nbr_atom = bond.GetOtherAtom(anomeric_carbon)
        if nbr_atom.GetAtomicNum() == 7:
            # Check if nitrogen is part of the nucleobase
            if mol.HasSubstructMatch(purine_pattern, useChirality=False, atomMap={nbr_atom.GetIdx(): nbr_atom.GetIdx()}):
                glycosidic_bond_found = True
                break
            if mol.HasSubstructMatch(pyrimidine_pattern, useChirality=False, atomMap={nbr_atom.GetIdx(): nbr_atom.GetIdx()}):
                glycosidic_bond_found = True
                break
    if not glycosidic_bond_found:
        return False, "No glycosidic bond between ribose and nucleobase found"

    return True, "Molecule is a ribonucleoside with D-ribose sugar and nucleobase connected via glycosidic bond"