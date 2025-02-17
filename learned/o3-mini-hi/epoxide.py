"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide
Definition: Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
An epoxide is defined here as a 3-membered ring that contains exactly one oxygen and two carbons,
all bonds in the ring must be single, and the carbons must not be directly conjugated (i.e. have
exocyclic double bonds) so that they represent a saturated epoxide motif.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_epoxide(smiles: str):
    """
    Determines whether the molecule (given as a SMILES string) can be classified as an epoxide.
    We define an epoxide as any 3-membered ring (cyclic) containing exactly one oxygen and two carbons,
    where all bonds in the ring are single and none of the carbon atoms is directly conjugated to an external
    double bond.

    Args:
        smiles (str): SMILES string for the molecule.

    Returns:
        bool: True if at least one valid epoxide ring is found, False otherwise.
        str: Detailed reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Retrieve ring information.
    ring_info = mol.GetRingInfo().AtomRings()  # each ring is a tuple of atom indices

    valid_epoxide_count = 0

    # Iterate over each ring.
    for ring in ring_info:
        if len(ring) != 3:
            continue  # Only consider 3-membered rings

        # Get atoms in the ring.
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Quick check: if any atom is aromatic, skip the ring.
        if any(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue

        # Count the types.
        oxygen_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 8)
        carbon_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 6)
        if oxygen_count != 1 or carbon_count != 2:
            continue

        # Check that all bonds in the ring are single bonds.
        # There are three bonds: between atom0-atom1, atom1-atom2, and atom2-atom0.
        bond01 = mol.GetBondBetweenAtoms(ring[0], ring[1])
        bond12 = mol.GetBondBetweenAtoms(ring[1], ring[2])
        bond20 = mol.GetBondBetweenAtoms(ring[2], ring[0])
        if not (bond01 and bond12 and bond20):
            continue
        if (bond01.GetBondType() != rdchem.BondType.SINGLE or
            bond12.GetBondType() != rdchem.BondType.SINGLE or
            bond20.GetBondType() != rdchem.BondType.SINGLE):
            continue

        # For each carbon in the epoxide ring, make sure that no bond going outside the ring is a double bond.
        ring_set = set(ring)
        skip_ring = False
        for atom in atoms_in_ring:
            if atom.GetAtomicNum() == 6:  # only check for carbons
                for bond in atom.GetBonds():
                    # if the neighbor is inside the ring, skip; else check bond type.
                    nei = bond.GetOtherAtom(atom)
                    if nei.GetIdx() in ring_set:
                        continue
                    # If an exocyclic double bond is found, then this ring is conjugated.
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        skip_ring = True
                        break
                if skip_ring:
                    break
        if skip_ring:
            continue

        # If all criteria are met, count this ring as a valid epoxide.
        valid_epoxide_count += 1

    if valid_epoxide_count > 0:
        return True, f"Found epoxide ring(s): {valid_epoxide_count} occurrence(s) of a 3-membered cyclic ether"
    else:
        return False, "No valid epoxide ring (3-membered cyclic ether with one oxygen and two carbons) found"


# (Optional test code)
if __name__ == "__main__":
    # A few test cases to see how classification works.
    test_smiles = [
        "CCCCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4; should be True
        "COc1cc2CN3CC[C@]4([C@H]5O[C@H]5CCC34)c2cc1OC",          # Zephyramine; expected NOT to be a valid isolated epoxide.
        "C1O[C@@H]1c1ccccc1",                                      # (R)-styrene oxide; should be True
    ]
    for s in test_smiles:
        result, reason = is_epoxide(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")