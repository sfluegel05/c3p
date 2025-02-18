"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Step 1: Find a ribose-like sugar (furanose ring with at least two hydroxyl groups)
    # More flexible pattern: O-C1-C-O with possible substitutions
    sugar_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O1")  # Basic furanose with multiple OHs
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        # Try alternative pattern for deoxyribose (missing one OH)
        sugar_pattern_alt = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H]([C@H](CO)O1)")
        sugar_matches = mol.GetSubstructMatches(sugar_pattern_alt)
        if not sugar_matches:
            return False, "No ribose-like sugar detected"

    # Step 2: Verify nitrogenous base attachment to C1' (first carbon in sugar pattern)
    base_attached = False
    for sugar_match in sugar_matches:
        c1_idx = sugar_match[0]  # First atom in sugar pattern is C1'
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        for neighbor in c1_atom.GetNeighbors():
            # Look for N-containing heterocycle (base)
            if neighbor.GetAtomicNum() == 7 and neighbor.IsInRing():
                base_attached = True
                break
        if base_attached:
            break
    if not base_attached:
        return False, "No nitrogenous base attached to sugar"

    # Step 3: Check for phosphate ester at 3' or 5' position
    phosphate_found = False
    # Phosphate connected via oxygen to sugar's 3' or 5' OH
    phosphate_pattern = Chem.MolFromSmarts("[O;!H0]-P(=O)([O-])[O]")

    # 5' position is typically the CH2OH group in the sugar pattern
    for sugar_match in sugar_matches:
        # 5' oxygen is the O in CO (last atom in sugar_match)
        o5_idx = sugar_match[-1]
        # 3' oxygen is the third hydroxyl in the pattern (index 2)
        o3_idx = sugar_match[2]

        # Check if either oxygen is connected to a phosphate
        for o_idx in [o5_idx, o3_idx]:
            atom = mol.GetAtomWithIdx(o_idx)
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    neighbor = bond.GetOtherAtomIdx(o_idx)
                    if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 15:  # Phosphorus
                        phosphate_found = True
                        break
            if phosphate_found:
                break
        if phosphate_found:
            break

    if not phosphate_found:
        return False, "No phosphate group at 3' or 5' position"

    # Step 4: Exclude CoA-like structures with thioester linkages
    if mol.HasSubstructMatch(Chem.MolFromSmarts("SCCNC(=O)")):
        return False, "Contains CoA-like thioester"

    return True, "Nucleoside with phosphate at 3' or 5' position"