"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:15727 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose/deoxyribose sugar with a phosphate group
    attached to the 5' carbon and a nitrogenous base (purine/pyrimidine) attached to the 1' carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group attached to 5' position (CH2OP)
    phosphate_pattern = Chem.MolFromSmarts("[CH2]OP(=O)")
    if phosphate_pattern is None:
        return False, "Invalid phosphate SMARTS pattern"
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group attached to 5' position"

    # Verify the CH2 is part of a 5-membered sugar ring with oxygen
    valid_phosphate = False
    for match in phosphate_matches:
        ch2_idx = match[0]
        ch2_atom = mol.GetAtomWithIdx(ch2_idx)
        # Check neighbors of CH2 for connection to sugar ring
        for neighbor in ch2_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.IsInRing():
                # Check if this neighbor is in a 5-membered ring with oxygen
                rings = mol.GetRingInfo().AtomRings()
                for ring in rings:
                    if neighbor.GetIdx() in ring and len(ring) == 5:
                        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                            valid_phosphate = True
                            break
                if valid_phosphate:
                    break
        if valid_phosphate:
            break
    if not valid_phosphate:
        return False, "Phosphate not attached to a 5-membered sugar ring"

    # Check for nitrogenous base (purine/pyrimidine) attached to 1' carbon of sugar
    # Look for glycosidic bond: N connected to sugar ring carbon
    glycosidic_found = False
    rings = mol.GetRingInfo().AtomRings()
    sugar_rings = [ring for ring in rings if len(ring) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)]
    for s_ring in sugar_rings:
        for atom_idx in s_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'N' and neighbor.IsInRing():
                        # Check if the nitrogen is part of a base ring with >=2 nitrogens
                        base_rings = mol.GetRingInfo().AtomRings()
                        for b_ring in base_rings:
                            if neighbor.GetIdx() in b_ring:
                                n_count = sum(1 for idx in b_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                                if n_count >= 2:
                                    glycosidic_found = True
                                    break
                        if glycosidic_found:
                            break
                if glycosidic_found:
                    break
        if glycosidic_found:
            break
    if not glycosidic_found:
        return False, "No nitrogenous base attached via glycosidic bond"

    return True, "Contains 5'-phosphate group attached to ribose/deoxyribose with nitrogenous base"