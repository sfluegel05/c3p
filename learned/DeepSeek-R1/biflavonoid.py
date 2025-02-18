"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, BondType

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer with two aryl-substituted benzopyran rings joined by a single bond or atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings in the molecule
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()

    benzopyran_subunits = []

    for ring in rings:
        # Check if the ring is a benzene ring (6 aromatic carbons)
        if len(ring) != 6:
            continue
        is_benzene = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_benzene:
            continue

        # Find fused rings that contain an oxygen
        fused_oxygen_rings = []
        for other_ring in rings:
            if other_ring == ring:
                continue
            shared = set(ring) & set(other_ring)
            if len(shared) >= 2:  # Fused rings share at least two atoms
                # Check if other_ring has at least one oxygen
                has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in other_ring)
                if has_oxygen:
                    fused_oxygen_rings.append(other_ring)
        if fused_oxygen_rings:
            # Take the first fused oxygen ring (assuming each benzene is part of one benzopyran)
            oxygen_ring = fused_oxygen_rings[0]
            subunit = set(ring).union(oxygen_ring)
            benzopyran_subunits.append(subunit)

    # Check for exactly two benzopyran subunits
    if len(benzopyran_subunits) != 2:
        return False, f"Found {len(benzopyran_subunits)} benzopyran subunits, need exactly 2"

    subunit1, subunit2 = benzopyran_subunits

    # Check if subunits are connected by a single bond or a shared atom
    connecting_bonds = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in subunit1 and a2 in subunit2) or (a1 in subunit2 and a2 in subunit1):
            connecting_bonds.append(bond)

    common_atoms = subunit1.intersection(subunit2)

    # Check connection criteria
    if len(common_atoms) == 1:
        return True, "Two benzopyran subunits connected via a shared atom"
    elif len(connecting_bonds) == 1 and connecting_bonds[0].GetBondType() == BondType.SINGLE:
        return True, "Two benzopyran subunits connected by a single bond"
    else:
        reasons = []
        if len(connecting_bonds) != 1:
            reasons.append(f"{len(connecting_bonds)} connecting bonds")
        else:
            if connecting_bonds[0].GetBondType() != BondType.SINGLE:
                reasons.append("non-single bond connection")
        if len(common_atoms) > 0:
            reasons.append(f"{len(common_atoms)} shared atoms")
        reason = "Subunits not properly connected: " + ", ".join(reasons) if reasons else "Subunits not properly connected"
        return False, reason