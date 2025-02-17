"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide
Definition: Any cyclic ether in which the oxygen atom forms part of a 3‐membered ring.
An epoxide has one oxygen and two carbons in a three‐membered ring, with all bonds being single bonds
and (typically) all atoms adopting sp³ hybridization.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_epoxide(smiles: str):
    """
    Determines whether the molecule (given as a SMILES string) can be classified as an epoxide.
    We define an epoxide as any 3-membered ring (cyclic) containing exactly one oxygen and two carbons,
    where all bonds in the ring are single bonds and all three atoms are sp³-hybridized.

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

    valid_epoxide_count = 0
    ring_info = mol.GetRingInfo().AtomRings()  # list of tuples of atom indices for each ring

    # Iterate over each ring and check if it is a 3-membered ring matching our criteria.
    for ring in ring_info:
        if len(ring) != 3:
            continue  # not a 3-membered ring

        # Get the atoms in the ring.
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Count how many oxygens and carbons are present.
        oxygen_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 8)
        carbon_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 6)
        # For a valid epoxide, we should have exactly one oxygen and two carbons.
        if oxygen_count != 1 or carbon_count != 2:
            continue

        # Check that the bonds connecting the 3 atoms are all single bonds.
        # In a 3-membered ring the bonds are between: atom0-atom1, atom1-atom2, and atom2-atom0.
        bond01 = mol.GetBondBetweenAtoms(ring[0], ring[1])
        bond12 = mol.GetBondBetweenAtoms(ring[1], ring[2])
        bond20 = mol.GetBondBetweenAtoms(ring[2], ring[0])
        if not (bond01 and bond12 and bond20):
            continue
        if (bond01.GetBondType() != rdchem.BondType.SINGLE or
            bond12.GetBondType() != rdchem.BondType.SINGLE or
            bond20.GetBondType() != rdchem.BondType.SINGLE):
            continue

        # Check that all atoms in the ring are sp³ hybridized.
        # (Epoxides are typically sp³; a deviation may indicate a different chemical environment.)
        if not all(atom.GetHybridization() == rdchem.HybridizationType.SP3 for atom in atoms_in_ring):
            continue

        # If all criteria are met, count this ring as a valid epoxide.
        valid_epoxide_count += 1

    if valid_epoxide_count > 0:
        return True, f"Found epoxide ring(s): {valid_epoxide_count} occurrence(s) of a 3-membered cyclic ether"
    else:
        return False, "No valid epoxide ring (3-membered cyclic ether with one oxygen and two carbons) found"

# (Optional test code)
if __name__ == "__main__":
    # Examples and counterexamples can be tried here.
    test_smiles = [
        "CCCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4 (should be True)
        "O=C1C2=C(O)C3=C(O)C=C(OC)C4=C3C(=C2C5(O)C6(C1(O)C(O)=C(C(=O)C6)C(=O)N)O5)CC74C(=CCCC7(C)C)C"  # example false positive
    ]
    for s in test_smiles:
        result, reason = is_epoxide(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")