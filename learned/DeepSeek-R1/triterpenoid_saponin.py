"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: CHEBI:????? triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import Mol, MolFromSmiles, MolFromSmarts
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import Atom

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin consists of a triterpenoid aglycone with one or more glycosyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for presence of glycosyl (sugar) groups
    ring_info = mol.GetRingInfo()
    sugar_found = False

    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue  # focus on common sugar ring sizes

        # Check if the ring contains an oxygen atom
        has_ring_oxygen = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if not has_ring_oxygen:
            continue

        # Count hydroxyl groups attached to ring carbons
        oh_count = 0
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetAtomicNum() != 6:
                continue  # only check carbons in the ring

            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    # Check if oxygen is part of a hydroxyl group (single bond)
                    if neighbor.GetDegree() == 1:
                        oh_count += 1

        if oh_count >= 2:
            sugar_found = True
            break

    if not sugar_found:
        return False, "No sugar group detected"

    # Check molecular weight (triterpenoid + sugar typically >500)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight ({mol_wt:.1f}) too low"

    # Check carbon count (aglycone ~30 + sugar(s))
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Only {c_count} carbons, insufficient for triterpenoid"

    # Check ring count (triterpenoid typically has >=4 rings)
    n_rings = len(ring_info.AtomRings())
    if n_rings < 4:
        return False, f"Only {n_rings} rings, insufficient for triterpenoid"

    return True, "Contains glycosyl group and meets triterpenoid criteria"