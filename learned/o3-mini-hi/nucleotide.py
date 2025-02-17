"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is defined as a nucleoside phosphate resulting from the condensation 
of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This function tests for:
   - Presence of at least one phosphorus atom (phosphate group)
   - A furanose sugar ring (5-membered ring with exactly 1 oxygen and 4 carbons)
   - A nucleobase-like aromatic heterocycle (an aromatic ring with at least 2 nitrogen atoms)
   - The phosphate group attached to the sugar ring.
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.

    A nucleotide should contain a nucleoside (a sugar ring attached to a nucleobase)
    with a phosphate group esterified at the 3' or 5' hydroxy position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a nucleotide, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorus atoms to confirm presence of a phosphate group.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found (missing phosphate group)"

    # Identify a furanose sugar ring: a 5-membered ring with 1 oxygen and 4 carbons.
    sugar_ring_found = False
    sugar_ring_indices = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            num_ox = 0
            num_c = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'O':
                    num_ox += 1
                elif atom.GetSymbol() == 'C':
                    num_c += 1
            if num_ox == 1 and num_c == 4:
                sugar_ring_found = True
                sugar_ring_indices = set(ring)
                break

    if not sugar_ring_found:
        return False, "No furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons) found"

    # Identify a nucleobase: search aromatic rings (other than the sugar ring) with at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        # Skip ring if it is identical to the sugar ring.
        if ring_set == sugar_ring_indices:
            continue
        # For fused systems the sugar ring may share atoms with the nucleobase;
        # so require at least one atom not belonging to the sugar ring.
        if ring_set.issubset(sugar_ring_indices):
            continue

        # Check if all atoms in the ring are aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen >= 2:
                nucleobase_found = True
                break

    if not nucleobase_found:
        return False, "No nucleobase-like aromatic heterocycle (ring with ≥2 nitrogens) found"

    # Confirm that at least one phosphate group is attached to the sugar ring.
    phosphate_attached = False
    for p_atom in p_atoms:
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetIdx() in sugar_ring_indices:
                phosphate_attached = True
                break
        if phosphate_attached:
            break

    if not phosphate_attached:
        return False, "Phosphate group is not attached to the sugar ring"

    return True, "Molecule contains a nucleoside portion (sugar and base) with a phosphate group attached"