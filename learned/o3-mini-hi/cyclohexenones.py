"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Any six-membered alicyclic ketone having one double bond in the ring (cyclohexenones).
Now uses a stricter enone check: the ketone group must be conjugated to the unique alkene in the ring.
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    
    A cyclohexenone is defined as any six-membered non-aromatic ring that has exactly one alkene (C=C)
    bond inside the ring AND a ketone functionality (C=O) attached to one of the alkene carbons (i.e. conjugated).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as cyclohexenone, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings as lists of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()

    # Loop over each ring in the molecule.
    for ring in ring_info:
        if len(ring) != 6:
            continue  # we only consider six-membered rings

        # Ensure the ring is non-aromatic
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Collect all bonds that lie completely inside the ring.
        ring_set = set(ring)
        alkene_bonds = []  # list of tuples (bond, atom_idx1, atom_idx2)
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                # Count double bonds only between two carbons.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    atom1 = mol.GetAtomWithIdx(a1)
                    atom2 = mol.GetAtomWithIdx(a2)
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        alkene_bonds.append((bond, a1, a2))

        # We expect exactly one alkene double bond within the ring.
        if len(alkene_bonds) != 1:
            continue

        # For the single alkene bond, check for an exocyclic ketone (C=O) conjugated to it.
        bond, a1_idx, a2_idx = alkene_bonds[0]
        def has_conjugated_ketone(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                # Require that the neighbor is not in the ring.
                if nbr.GetIdx() in ring_set:
                    continue
                # Look for an oxygen double-bonded to the atom.
                if nbr.GetAtomicNum() == 8:
                    ex_bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if ex_bond is not None and ex_bond.GetBondType() == Chem.BondType.DOUBLE:
                        return True
            return False

        if has_conjugated_ketone(a1_idx) or has_conjugated_ketone(a2_idx):
            return True, ("Found a non-aromatic six-membered ring with exactly one alkene bond " 
                          "conjugated to a ketone carbonyl (cyclohexenone motif).")
        else:
            continue  # the ring did not show the conjugated ketone requirement

    return False, "No six-membered ring with a conjugated alkene and ketone (cyclohexenone) was found."

# Example test calls (uncomment to run)
# test_smiles = "O[C@H]1[C@H](O)C(=CC(=O)C1N)CO"  # 2-Aminovalienone example
# result, reason = is_cyclohexenones(test_smiles)
# print(result, reason)