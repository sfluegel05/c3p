"""
Classifies: CHEBI:17002 cholesteryl ester
"""
#!/usr/bin/env python
"""
Classifies: Cholesteryl esters

Definition: A sterol ester obtained by formal condensation of the carboxy group
of any carboxylic acid with the 3-hydroxy group of cholesterol.
In cholesteryl esters the free 3-hydroxy in cholesterol is esterified,
which means that the oxygen of the diminished OH group is now connected to 
an acyl group (via a carbonyl).
This implementation searches for an ester bond with the oxygen of the bond 
attached to a cyclic (preferably multi‐ring) system. We also require that the 
molecule contains at least four rings (the steroid nucleus is a tetracyclic fused ring system).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The method is twofold:
      1. Find an ester bond: an oxygen (O_est) bonded to a carbonyl carbon (C=O)
         (i.e. a single O–C bond where that carbon is double bonded to an oxygen).
      2. Verify that the O_est is connected (via another bond) to a cyclic (ring) system,
         and that the molecule contains at least four rings (indicating a steroid nucleus).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is identified as a cholesteryl ester, False otherwise.
        str : Explanation of the decision.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if there are at least four rings in the overall molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 4:
        return False, f"Found only {len(rings)} rings; cholesteryl esters require a tetracyclic (4-ring) steroid nucleus"
    
    # Iterate over bonds looking for an ester linkage.
    for bond in mol.GetBonds():
        # We are looking for a single bond between an oxygen and a carbon.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify which atom is oxygen and which is carbon.
        if a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6:
            o_atom = a1
            c_atom = a2
        elif a2.GetAtomicNum() == 8 and a1.GetAtomicNum() == 6:
            o_atom = a2
            c_atom = a1
        else:
            continue
        
        # Check that the carbon (c_atom) is part of a carbonyl group (C=O),
        # i.e. it is double-bonded to an oxygen (other than the one in the ester bond).
        has_carbonyl = False
        for nb in c_atom.GetNeighbors():
            if nb.GetIdx() == o_atom.GetIdx():
                continue  # skip the oxygen of the ester bond
            # Check if the neighbor is oxygen and the bond is a double bond.
            bond_to_nb = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nb.GetIdx())
            if nb.GetAtomicNum() == 8 and bond_to_nb.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if not has_carbonyl:
            continue  # not an ester bond
        
        # At this point we have an ester-like O-C(=O) bond.
        # Now check that the oxygen (which came from the alcohol part) 
        # is attached (by another bond) to a ring system.
        attached_to_ring = False
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue  # skip the carbonyl carbon
            if nbr.IsInRing():
                attached_to_ring = True
                break
        if not attached_to_ring:
            continue  # the oxygen is not attached to a cyclic (steroid like) moiety
        
        # If we reach here, we found an ester bond whose oxygen is attached to a ring,
        # and the overall molecule has at least four rings (indicative of a steroid nucleus).
        return True, ("Contains an ester linkage (O–C(=O)) with the oxygen attached to a cyclic system "
                      "and the molecule contains at least 4 rings (a steroid nucleus), "
                      "consistent with a cholesteryl ester")
    
    return False, "No ester bond attached to a polycyclic (steroid-like) nucleus found; not a cholesteryl ester"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)[C@H](C)CCCC(C)C"  # cholesteryl linoleate
    result, reason = is_cholesteryl_ester(test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)