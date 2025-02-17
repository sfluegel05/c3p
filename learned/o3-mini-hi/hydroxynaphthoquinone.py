"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety 
            (a fused bicyclic ring system that is “naphthalene‐like” and
             contains two carbonyl (C=O) groups attached to two of the 10 carbons)
            is substituted by at least one hydroxy group directly attached 
            to one of the ring carbons.

This implementation searches for fused pairs of aromatic six‐membered rings 
that share at least two atoms and with a union of 10 atoms (like naphthalene).
It then requires that exactly two of these 10 atoms be carbonyl carbons (C=O)
and that at least one atom of the candidate core has an external –OH substituent.

Note: This heuristic may miss molecules where the naphthoquinone “core” is 
extended (i.e. fused to additional rings) or altered, but it improves on false positive
classification by demanding an exact match to a 10–atom naphthalene system.
"""

from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone.
    
    Algorithm:
      1. Parse the SMILES string.
      2. Get all ring systems from the molecule.
      3. Look for pairs of aromatic rings (each of size 6) that share at least
         two atoms and whose union is exactly 10 atoms (i.e. a naphthalene-like core).
      4. In each candidate fused system, count how many atoms are involved in a carbonyl group.
         A carbonyl is defined here as an aromatic carbon atom that has at least one double bond 
         to an oxygen (atomic number 8).
      5. For each candidate core with exactly 2 carbonyl substituents, check if at least one atom
         in the core has an external substituent that is –OH (an oxygen atom bound via a single bond, 
         carrying at least one hydrogen).
      6. If a candidate passes the above tests, we return True plus an explanation.
      
    Args:
         smiles (str): SMILES string representation of the molecule.
    
    Returns:
         (bool, str): Tuple with a True/False classification and a reason.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices

    candidate_found = False
    reason_msg = ""
    
    # Look for pairs of rings that meet the naphthalene core criteria.
    # We require each ring to be of size 6 and aromatic.
    for i in range(len(atom_rings)):
        ring1 = atom_rings[i]
        if len(ring1) != 6:
            continue
        # check ring1 atoms are aromatic carbons
        if any(not mol.GetAtomWithIdx(idx).GetIsAromatic() or mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring1):
            continue
        for j in range(i+1, len(atom_rings)):
            ring2 = atom_rings[j]
            if len(ring2) != 6:
                continue
            if any(not mol.GetAtomWithIdx(idx).GetIsAromatic() or mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring2):
                continue
            
            # Check for fusion: these two rings should share at least 2 atoms.
            common = set(ring1).intersection(ring2)
            if len(common) < 2:
                continue
            
            # The candidate naphthalene-like core is the union of both rings.
            core_atoms = set(ring1).union(ring2)
            if len(core_atoms) != 10:
                # If the union is larger than 10 then other rings are attached 
                # or the system is extended.
                continue
            
            # Now check that all atoms in the candidate core are aromatic carbons.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in core_atoms):
                continue
            
            # In a proper naphthoquinone, two of these carbons should be converted to carbonyls.
            # Count carbonyl groups in the core.
            carbonyl_count = 0
            for idx in core_atoms:
                atom = mol.GetAtomWithIdx(idx)
                # Look at bonds: a double bond to an oxygen (atomic number 8)
                for bond in atom.GetBonds():
                    # Ensure the bond is double and the other atom is oxygen.
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        nbr = bond.GetOtherAtom(atom)
                        if nbr.GetAtomicNum() == 8:
                            carbonyl_count += 1
                            break  # count each carbon only once
            
            if carbonyl_count != 2:
                # Not the typical naphthoquinone pattern.
                continue
            
            # Now check for hydroxy substituents attached directly to the core.
            # For each atom in the candidate, look for neighbors (outside core) that 
            # are single-bonded oxygen with at least one hydrogen.
            hydroxy_count = 0
            for idx in core_atoms:
                atom = mol.GetAtomWithIdx(idx)
                for bond in atom.GetBonds():
                    # only consider single bonds for substituents 
                    if bond.GetBondType() != Chem.BondType.SINGLE:
                        continue
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetIdx() in core_atoms:
                        continue
                    # Check that the neighbor is oxygen and has at least one H.
                    if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                        hydroxy_count += 1
            if hydroxy_count >= 1:
                candidate_found = True
                reason_msg = (f"Found a naphthoquinone core (10 atoms with 2 carbonyls) with "
                              f"{hydroxy_count} hydroxy substituent(s) attached")
                return True, reason_msg
    if not candidate_found:
        # We differentiate between “no naphthoquinone core found” versus “core found but no OH substituents”
        # For simplicity, if no candidate core meeting the criteria was found, we return a generic message.
        return False, "No valid 10–atom naphthoquinone core with 2 carbonyl groups and a hydroxy substituent found"
    
    return False, "Could not classify the molecule as hydroxynaphthoquinone"

# (Optional) Uncomment the code below to test the function with sample SMILES.
# if __name__ == "__main__":
#     test_smiles = [
#         # True positive examples:
#         ("COC1=C(C)C(=O)c2c(O)cc(OC\\C=C(/C)CCC=C(C)C)cc2C1=O", "7-O-geranyl-2-O,3-dimethylflaviolin"),
#         ("CC(=O)OC(CC=C(C)C)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "Acetylshikonin"),
#         ("Oc1ccc(O)c2C(=O)C=CC(=O)c12", "naphthazarin"),
#         ("Oc1cccc2C(=O)C=CC(=O)c12", "juglone"),
#         ("OC1=CC(=O)c2ccccc2C1=O", "lawsone"),
#         ("Cc1cc(O)c2C(=O)C=CC(=O)c2c1", "Ramentaceone"),
#         ("C1=CC=C(C2=C1C(C=C(C2=O)O)=O)O", "2,8-dihydroxy-1,4-naphthoquinone"),
#         # False negative example:
#         ("[C@H](C)([C@@H]([C@@H]([C@H](\\C=C\\O[C@]1(OC=2C(C1=O)=C3C(C(C(=C(C3=O)/C=N/N4CCN(CC4)C)[O-])=O)=C(C2C)[O-])C)OC)C)OC(=O)C)[C@H](O)[C@@H]([C@@H](O)[C@@H](C)/C=C/C=C(/C)C(N)=O)C",
#          "rifampicin para-naphthoquinone carboxamide(2-)"),
#     ]
#     for smi, name in test_smiles:
#         res, msg = is_hydroxynaphthoquinone(smi)
#         print(f"NAME: {name}\nSMILES: {smi}\nResult: {res}\nReason: {msg}\n")