"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N‐acylsphingosine (parent compounds of the ceramide family)
Definition: composed of a sphingosine backbone (i.e. an acyclic chain featuring a secondary amine 
with two hydroxyl‐bearing carbons) with an unspecified fatty acyl group attached to the nitrogen.
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is defined as a sphingosine backbone (an open-chain structure containing
    a nitrogen attached to at least two hydroxylated carbons) in which the nitrogen is acylated 
    (attached to a fatty acyl group via an amide bond).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule matches the N-acylsphingosine criteria, False otherwise.
        str: A message indicating the reason for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the sphingosine backbone.
    # We expect a nitrogen attached to two carbons that bear hydroxyl groups.
    # The chiral indicators ([C@@H] and [C@H]) are common in sphingosine examples.
    sphingosine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@H](O)")
    if sphingosine_pattern is None:
        return False, "Error in sphingosine SMARTS definition"

    # Find matches for the sphingosine backbone in the molecule.
    backbone_matches = mol.GetSubstructMatches(sphingosine_pattern)
    if not backbone_matches:
        return False, "Sphingosine backbone not found"

    # For each match, check two additional criteria:
    # 1. The matched atoms (nitrogen and the two carbons in the backbone pattern) should all be acyclic.
    #    This distinguishes sphingosine from similar patterns embedded in sugars or ring systems.
    # 2. The nitrogen (first atom in the match) should be acylated.
    #    That is, aside from the backbone connection, it should be bonded to a carbon
    #    that is part of a carbonyl group (i.e. has a double bond to an oxygen).
    for match in backbone_matches:
        # Check that none of the backbone atoms are in a ring
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            continue  # Skip this match if any backbone atom is part of a ring

        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Get the indices of atoms connected to the nitrogen.
        neighbor_indices = [nbr.GetIdx() for nbr in n_atom.GetNeighbors()]
        acyl_found = False
        for nbr_idx in neighbor_indices:
            # Skip the backbone neighbor (match[1] is the carbon from the backbone)
            if nbr_idx == match[1]:
                continue
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            # We expect the acyl group to be attached via an amide bond.
            # So the neighboring atom should be a carbon bonded to an oxygen by a double bond.
            if nbr_atom.GetAtomicNum() == 6:  # carbon
                for bond in nbr_atom.GetBonds():
                    # Look for a bond from this carbon to an oxygen that is double-bonded.
                    other_atom = bond.GetOtherAtom(nbr_atom)
                    if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        acyl_found = True
                        break
            if acyl_found:
                break
        if acyl_found:
            return True, "Molecule contains a sphingosine backbone (acyclic) with an N-linked acyl (fatty acid) group"
    
    return False, "No acylated sphingosine backbone meeting all criteria found"