"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base attached to a deoxyribose sugar via a β-N1-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and small fragments
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Exclude molecules with phosphate groups (nucleotides)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Molecule is a nucleotide (contains phosphate group)"

    # Define pyrimidine base SMARTS (allows for substitutions)
    pyrimidine_smarts = "n1c([aH0])[nH0]c([aH0])c([aH0])c1[aH0]"
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)
    if mol.HasSubstructMatch(pyrimidine_pattern):
        pyrimidine_match = mol.GetSubstructMatch(pyrimidine_pattern)
    else:
        return False, "No pyrimidine base found"

    # Define deoxyribose sugar SMARTS (five-membered ring with oxygen and without 2' OH)
    deoxyribose_smarts = "[C@@H]1O[C@@H]([C@H](C1)O)CO"
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_smarts)
    if mol.HasSubstructMatch(deoxyribose_pattern):
        sugar_match = mol.GetSubstructMatch(deoxyribose_pattern)
    else:
        return False, "No deoxyribose sugar found"

    # Check for glycosidic bond between N1 of pyrimidine and C1' of sugar
    n1_atom_idx = pyrimidine_match[0]  # Assuming N1 is the first atom in the pattern
    c1_prime_idx = sugar_match[0]      # Assuming C1' is the first atom in the sugar pattern
    bond = mol.GetBondBetweenAtoms(n1_atom_idx, c1_prime_idx)
    if bond is None:
        return False, "No glycosidic bond between pyrimidine base N1 and sugar C1'"

    # Ensure the molecule is not overly complex (no additional rings)
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    if num_rings > 2:
        return False, f"Molecule has {num_rings} rings; expected 2 (pyrimidine and sugar)"

    # Check for additional substituents on the sugar (excluding necessary hydroxyls)
    sugar_atoms = set(sugar_match)
    for atom_idx in sugar_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in sugar_atoms and neighbor.GetAtomicNum() != 8:
                    return False, "Sugar has additional substituents"

    # Check for additional substituents on the pyrimidine base
    base_atoms = set(pyrimidine_match)
    for atom_idx in base_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() in (6, 7):  # Carbon or nitrogen atom
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in base_atoms and neighbor.GetAtomicNum() != 1:
                    return False, "Pyrimidine base has additional substituents"

    return True, "Contains pyrimidine base attached to deoxyribose via N-glycosidic bond"