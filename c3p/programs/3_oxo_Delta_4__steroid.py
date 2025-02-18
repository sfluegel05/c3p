"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-Delta(4) steroid 
Definition: A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position,
i.e. a steroid (with a fused tetracyclic nucleus) that contains an α,β-unsaturated ketone 
(enone) in one of its six-membered rings.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    The molecule must have a fused steroid nucleus (three six-membered rings and one five-membered ring,
    with a sufficient number of ring atoms) and, in one of the six-membered rings, an enone motif
    is identified. The enone is defined here as a carbon that bears a carbonyl (C=O) group and is 
    connected by a single bond to a carbon that is part of a C=C double bond with a neighboring carbon (all within the same ring).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------------------------
    # 1. Verify steroid nucleus structure.
    # ---------------------------
    # Steroids typically have a fused tetracyclic nucleus consisting of three six-membered rings
    # and one five-membered ring. We use ring information as a heuristic.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Count six- and five-membered rings.
    count6 = sum(1 for r in rings if len(r) == 6)
    count5 = sum(1 for r in rings if len(r) == 5)
    if count6 < 3 or count5 < 1:
        return False, "Steroid nucleus not found (insufficient six-membered and/or five-membered rings)"
    
    # Build a set of all atoms that are in rings of size 5 or 6.
    steroid_atoms = set()
    for r in rings:
        if len(r) in [5,6]:
            steroid_atoms.update(r)
    if len(steroid_atoms) < 15:
        return False, "Steroid nucleus not found (fused ring atom count too low)"
    
    # ---------------------------
    # 2. Detect enone functionality within one six-membered ring of the nucleus.
    # ---------------------------
    # Instead of a fixed SMARTS pattern, we search each six-membered ring for an enone motif.
    # The desired pattern is:
    #   (a) A carbon atom (C1) that is part of a six-membered ring and has a double bond to an oxygen (C=O).
    #   (b) C1 must be connected (via a single bond) to a beta carbon (C2) in the ring.
    #   (c) This beta carbon C2 must participate in a double bond (C=C) with a gamma carbon (C3) in the same ring.
    rings6 = [set(r) for r in rings if len(r) == 6]
    
    enone_found = False
    for ring in rings6:
        # For each atom index in the six-membered ring, check for the pattern.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look for carbonyl carbon: must be carbon and have a double bond to oxygen.
            if atom.GetAtomicNum() != 6:
                continue
            has_carbonyl = False
            # List to hold the oxygen neighbor indices involved in C=O.
            carbonyl_oxygens = []
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        has_carbonyl = True
                        carbonyl_oxygens.append(nbr.GetIdx())
            if not has_carbonyl:
                continue
            
            # Now look for a beta carbon: a neighbor (other than the oxygen) that is in the same ring.
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetIdx() not in ring:
                    continue  # Only consider atoms in the same ring.
                if nbr.GetAtomicNum() != 6:
                    continue  # Should be a carbon.
                if nbr.GetIdx() in carbonyl_oxygens:
                    continue  # Skip the oxygen atom.
                
                # Check if this beta carbon participates in a C=C double bond
                # with another carbon (gamma carbon) also in the ring and not the original carbonyl carbon.
                for bond2 in nbr.GetBonds():
                    if bond2.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                        continue
                    gamma = bond2.GetOtherAtom(nbr)
                    if gamma.GetIdx() == atom.GetIdx():
                        continue  # Skip bond back to the carbonyl carbon.
                    if gamma.GetIdx() not in ring:
                        continue  # Must be in the same ring.
                    if gamma.GetAtomicNum() != 6:
                        continue
                    # Found a candidate enone: C1 (carbonyl) - C2 (beta, single bond) and C2=C3 (double bond).
                    enone_found = True
                    break
                if enone_found:
                    break
            if enone_found:
                break
        if enone_found:
            break
            
    if not enone_found:
        return False, "Enone motif not found in a six-membered ring within the steroid nucleus"
    
    return True, "Contains a fused steroid nucleus with a 3-oxo/Δ(4) enone motif in a six-membered ring."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: (20S)-20-hydroxypregn-4-en-3-one
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)O"
    result, reason = is_3_oxo_Delta_4__steroid(test_smiles)
    print(result, reason)