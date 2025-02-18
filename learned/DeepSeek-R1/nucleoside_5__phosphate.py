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

    # Check for phosphate group attached to 5' position of sugar
    # Pattern: Phosphate connected via ester to CH2 group that is part of a 5-membered ring with oxygen
    # SMARTS: [PX4](=O)(O)(O)OC[C@H]1O (adjusting for possible stereochemistry)
    phosphate_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[PX4](=O)(O)(O)OC[C@H]1O"))
    if not phosphate_matches:
        return False, "No phosphate group attached to 5' position of sugar"

    # Verify the connected carbon (4') is in a 5-membered ring with oxygen
    valid_phosphate = False
    for match in phosphate_matches:
        # The fourth atom in the match is the 4' carbon (C in OC[C@H]1O)
        c4_idx = match[3]
        c4_atom = mol.GetAtomWithIdx(c4_idx)
        # Check if this carbon is in a 5-membered ring with oxygen
        rings = mol.GetRingInfo().AtomRings()
        for ring in rings:
            if c4_idx in ring and len(ring) == 5:
                # Check for oxygen in the ring
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                    valid_phosphate = True
                    break
        if valid_phosphate:
            break
    if not valid_phosphate:
        return False, "Phosphate not attached to a 5-membered sugar ring with oxygen"

    # Check for nitrogenous base (purine/pyrimidine) attached to 1' carbon of sugar
    # Look for glycosidic bond: N-C1O (base N connected to sugar's 1' carbon)
    glycosidic_bond = Chem.MolFromSmarts("[NX3][C@H]1O[C@H]")
    if not mol.HasSubstructMatch(glycosidic_bond):
        return False, "No glycosidic bond (base-sugar linkage) found"

    # Verify the base has at least two nitrogens in a ring
    base_ring_nitrogens = 0
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        ring_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if ring_nitrogens >= 2:
            # Check if any nitrogen in the ring is connected to the sugar
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C' and neighbor.IsInRing():
                            # Check if neighbor is part of the sugar ring
                            neighbor_rings = mol.GetRingInfo().AtomRings()
                            for n_ring in neighbor_rings:
                                if neighbor.GetIdx() in n_ring and any(mol.GetAtomWithIdx(r_idx).GetAtomicNum() == 8 for r_idx in n_ring):
                                    base_ring_nitrogens += ring_nitrogens
                                    break
                            if base_ring_nitrogens >= 2:
                                break
                    if base_ring_nitrogens >= 2:
                        break
            if base_ring_nitrogens >= 2:
                break
    if base_ring_nitrogens < 2:
        return False, "Insufficient nitrogen atoms in base rings"

    return True, "Contains 5'-phosphate group attached to ribose/deoxyribose with nitrogenous base"