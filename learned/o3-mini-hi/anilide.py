"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide – Any aromatic amide obtained by acylation of aniline.

An anilide is defined as an aromatic amide that arises when an aniline (PhNH2) 
is acylated. In such molecules the amide nitrogen (the N–C(=O) part) originates from 
aniline so that one substituent on N is (or stems from) a benzene ring. 
This implementation scans for N–C(=O) bonds and then verifies that the nitrogen has,
aside from the acyl-carbon, at least one neighbor that is part of a benzene ring.
"""

from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide (obtained by acylation of aniline)
    based on its SMILES string.

    This function looks for an amide bond where the nitrogen is acylated (bonded to a carbonyl carbon)
    and also is connected to at least one benzene-like carbon (i.e. a carbon that belongs to a six-membered 
    aromatic ring composed solely of carbon atoms). We also require that the N atom carries at least one hydrogen,
    as expected for an acylation product of aniline (PhNH2).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an anilide, False otherwise.
        str: Reason for classification.
    """

    # Helper: returns True if the atom (which must be a carbon) is in a benzene ring:
    def is_benzene_carbon(atom, mol):
        if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
            return False
        ring_info = mol.GetRingInfo().AtomRings()
        for ring in ring_info:
            if atom.GetIdx() in ring and len(ring) == 6:
                # Ensure that every atom in this ring is a carbon and aromatic.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    return True
        return False

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        # Look for bonds between a nitrogen (candidate) and a carbon
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify which atom is nitrogen and which is carbon
        if a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            N_atom = a1
            C_atom = a2
        elif a2.GetAtomicNum() == 7 and a1.GetAtomicNum() == 6:
            N_atom = a2
            C_atom = a1
        else:
            continue

        # We want an amide bond: the carbon should be a carbonyl carbon (i.e. have a double bond to oxygen)
        carbonyl_found = False
        for nbond in C_atom.GetBonds():
            if nbond.GetBondType() == Chem.BondType.DOUBLE:
                other = nbond.GetOtherAtom(C_atom)
                if other.GetAtomicNum() == 8:
                    carbonyl_found = True
                    break
        if not carbonyl_found:
            continue

        # We now have an N–C(=O) bond. Next, check the environment of the nitrogen.

        # Check that the nitrogen has (at least) one neighbor (besides the carbonyl carbon) that belongs to a benzene ring.
        benzene_neighbor_found = False
        for neighbor in N_atom.GetNeighbors():
            # Skip the acyl carbon that provided the carbonyl
            if neighbor.GetIdx() == C_atom.GetIdx():
                continue
            if is_benzene_carbon(neighbor, mol):
                benzene_neighbor_found = True
                break

        if not benzene_neighbor_found:
            continue

        # Optional: require that the N atom has at least one hydrogen.
        # In an acylated aniline, the nitrogen originally comes from PhNH2.
        num_explicit_H = N_atom.GetTotalNumHs()
        if num_explicit_H < 1:
            # If no hydrogen, then this might be a tertiary amide: not classically from aniline.
            continue

        # If we reached here then:
        # - There is an N–C(=O) bond
        # - The nitrogen has at least one neighbor (besides the acyl group) that is a benzene carbon.
        # - The nitrogen still bears at least one hydrogen.
        return True, "Found anilide: amide N (with at least one H) attached to a benzene carbon and to a carbonyl group."

    return False, "No valid anilide substructure found"

# Test the function on one of the true positive examples:
if __name__ == "__main__":
    test_smiles = "CC(=O)NC1=CC=C(C=C1)C(=S)NCC2=CC=CO2"
    result, reason = is_anilide(test_smiles)
    print("Is anilide:", result)
    print("Reason:", reason)