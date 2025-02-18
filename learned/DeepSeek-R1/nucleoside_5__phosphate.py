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

    # Find all 5-membered sugar rings (ribose/deoxyribose)
    sugar_rings = []
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) != 5:
            continue
        # Check for 4 carbons and 1 oxygen (ribose) or 4 carbons and 0 oxygen (deoxyribose)
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        c_count = sum(1 for a in atoms if a.GetAtomicNum() == 6)
        o_count = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        if (c_count == 4 and o_count == 1) or (c_count == 5 and o_count == 0 and any(a.GetHybridization() == Chem.HybridizationType.SP3 for a in atoms)):
            sugar_rings.append(ring)

    if not sugar_rings:
        return False, "No ribose/deoxyribose ring found"

    # Check for phosphate group attached to 5' position (CH2OP connected to C4' of sugar)
    phosphate_found = False
    for s_ring in sugar_rings:
        # Find C4' (connected to CH2OP)
        for atom_idx in s_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'C':
                continue
            # Check if this atom (C4') is connected to a CH2OP group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() >= 2:
                    # Check for O-P attached to this CH2
                    for phosphate in neighbor.GetNeighbors():
                        if phosphate.GetSymbol() == 'O':
                            for p_atom in phosphate.GetNeighbors():
                                if p_atom.GetSymbol() == 'P':
                                    phosphate_found = True
                                    break
                            if phosphate_found:
                                break
                    if phosphate_found:
                        break
            if phosphate_found:
                break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group attached to 5' position of sugar"

    # Check for glycosidic bond between C1' (anomeric carbon) and nitrogenous base
    glycosidic_found = False
    for s_ring in sugar_rings:
        for atom_idx in s_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'C':
                continue
            # Check if this is the anomeric C (connected to O in ring and to N)
            is_anomeric = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.IsInRing():
                    is_anomeric = True
                    break
            if not is_anomeric:
                continue
            # Check if connected to a nitrogen in a heterocycle (base)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and not neighbor.IsInRingSize(5) and not neighbor.IsInRingSize(6):
                    # Check if the N is part of a purine/pyrimidine
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