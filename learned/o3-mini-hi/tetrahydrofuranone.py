"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: Tetrahydrofuranone
Definition: Any oxolane (i.e. a 5‐membered ring with one oxygen) having an oxo- substituent (a C=O group) 
at any position on the tetrahydrofuran ring.
"""

from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    A tetrahydrofuranone is defined as any molecule that contains a 5-membered ring (oxolane)
    with exactly one oxygen and four carbon atoms, where at least one of the carbon atoms in that ring 
    has a double-bonded oxygen (an oxo group) attached (which can be either part of a lactone ring or a substituent).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get information about rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples, each is a ring (list of atom indices)
    
    # Iterate over all rings
    for ring in rings:
        if len(ring) == 5:
            # Count the number of oxygen atoms in the ring:
            oxygen_count = 0
            carbon_indices = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxygen_count += 1
                elif atom.GetAtomicNum() == 6:
                    carbon_indices.append(idx)
            # Check if the ring qualifies as an oxolane (5-membered ring with exactly one oxygen)
            if oxygen_count != 1:
                continue  # not the tetrahydrofuran ring we are after
            
            # Now check if any carbon in the ring has an oxo group (i.e. a double-bonded oxygen)
            for c_idx in carbon_indices:
                atom = mol.GetAtomWithIdx(c_idx)
                for bond in atom.GetBonds():
                    # Check that the bond is a double bond and the neighbor is oxygen.
                    # We also want to ensure that the oxygen is not part of the ring itself (exocyclic oxo).
                    if bond.GetBondTypeAsDouble() == 2.0:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8:
                            # Optionally, we can check if the neighbor is outside the ring.
                            # However, in lactones the carbonyl oxygen is exocyclic relative to the ring framework.
                            # If needed, one could enforce: if neighbor.GetIdx() not in ring:
                            return True, f"Found tetrahydrofuranone: ring (atom indices {ring}) has an oxo group at carbon index {c_idx}."
            
    return False, "No tetrahydrofuranone substructure found."
    
# Examples (if you want to test):
if __name__ == "__main__":
    examples = [
        ("N-butyryl-L-homoserine lactone", "CCCC(=O)N[C@@H]1CCOC1=O"),
        ("gamma-decalactone", "C1(OC(=O)CC1)CCCCCC"),
        ("L-homoserine lactone", "N[C@H]1CCOC1=O"),
        ("gamma-caprolactone", "CCC1CCC(=O)O1"),
        ("Non-tetrahydrofuranone", "CC(C)CC")  # simple alkane fragment
    ]
    
    for name, smi in examples:
        res, reason = is_tetrahydrofuranone(smi)
        print(f"{name}: {res}, Reason: {reason}")