"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Any six-membered alicyclic ketone having one double bond in the ring (cyclohexenones)
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any non-aromatic six-membered ring that has one alkene double bond 
    within the ring and bears a ketone functionality (i.e. a carbonyl group on a ring carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (list of tuples of atom indices for each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    # Loop over all rings to search for a six-membered candidate
    for ring in ring_info:
        if len(ring) != 6:
            continue  # Only consider six-membered rings

        # Check that none of the ring atoms are aromatic (alicyclic requirement).
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Count alkene double bonds inside the ring (only count C=C bonds).
        alkene_double_count = 0
        # To avoid counting bonds twice, iterate over bonds that have both endpoints in the ring.
        ring_set = set(ring)
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            # Only consider bonds that connect two atoms in the ring.
            if begin_idx in ring_set and end_idx in ring_set:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Only count it if both atoms are carbons.
                    atom1 = bond.GetBeginAtom()
                    atom2 = bond.GetEndAtom()
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        alkene_double_count += 1

        # We expect exactly one alkene double bond inside the ring.
        if alkene_double_count != 1:
            continue
        
        # Check for the presence of a ketone on the ring: i.e. a ring carbon with a double bond to O.
        ketone_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check only carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors not necessarily in the ring. If there is an oxygen attached via a double bond,
            # consider it as a ketone group.
            for nbr in atom.GetNeighbors():
                # Skip if the neighbor is in the ring because we are looking for an exocyclic oxygen.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        ketone_found = True
                        break
            if ketone_found:
                break

        if not ketone_found:
            continue

        # If a six-membered ring meets both conditions, then the molecule is classified as a cyclohexenone.
        return True, "Found a non-aromatic six-membered ring with exactly one alkene bond and a ketone carbonyl attached."

    return False, "No six-membered non-aromatic ring with one alkene and a ketone group was found."