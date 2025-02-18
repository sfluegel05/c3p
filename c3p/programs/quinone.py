"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone 
Definition: Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones,
derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups
(with any necessary rearrangement of double bonds). (Polycyclic and heterocyclic analogues are included.)
"""

from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is assumed to be a compound with an aromatic ring that has at least two carbonyl (C=O) groups
    directly attached to atoms that are part of the ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a quinone, otherwise False.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information using GetRingInfo, which returns a tuple of atom index lists for each ring.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop through each ring
    for ring in rings: 
        # First check: Does the ring have at least 4 atoms? (Most aromatic rings are 5- or 6-membered)
        if len(ring) < 4:
            continue

        # Check if the ring is aromatic: all atoms in the ring should be flagged as aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Count how many atoms in the ring have a double bond to an oxygen (i.e. are part of a carbonyl group).
        carbonyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check only carbon atoms (atomic number 6)
            if atom.GetAtomicNum() != 6:
                continue
            # Check bonds from this atom to see if any are double bonds to oxygen.
            for bond in atom.GetBonds():
                # Ensure that the bond is a double bond and that the neighboring atom is oxygen (atomic number 8)
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        # Sometimes the oxygen might be counted more than once if the element appears in 
                        # multiple bonds. Once we found a double-bonded oxygen, count it and break.
                        carbonyl_count += 1
                        break

        # If we find at least two such carbonyl groups in one aromatic ring, we identify a quinone pattern.
        if carbonyl_count >= 2:
            return True, f"Found aromatic ring (atoms: {ring}) with {carbonyl_count} carbonyl groups."

    # If no qualifying ring was found, then the molecule is not a quinone.
    return False, "No aromatic ring with at least two conjugated carbonyl groups was found."

# Example calls (you can uncomment these lines to test the function)
# print(is_quinone("Oc1ccc2C(=O)c3c(O)ccc(O)c3C(=O)c2c1O"))  # quinalizarin
# print(is_quinone("O=C1C(OC)=CC(=O)C=C1"))  # a pattern similar to a benzoquinone