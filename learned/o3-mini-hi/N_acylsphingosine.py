"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N‐acylsphingosine (parent compounds of the ceramide family)
Definition: composed of a sphingosine backbone having an unspecified fatty acyl group attached to the nitrogen.
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine features a sphingosine backbone (a nitrogen linked to two hydroxyl-bearing carbons) 
    with an acyl (fatty acid) group attached to the nitrogen via an amide bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches the N-acylsphingosine criteria, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for the sphingosine backbone.
    # This pattern looks for a nitrogen attached to two carbons bearing hydroxyl groups,
    # which is characteristic of the sphingosine core.
    # Note: The stereochemistry markers [C@@H] and [C@H] are common in sphingosine examples.
    sphingosine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@H](O)")
    if sphingosine_pattern is None:
        return False, "Error in sphingosine SMARTS definition"

    # Find matches for the sphingosine backbone in the molecule.
    backbone_matches = mol.GetSubstructMatches(sphingosine_pattern)
    if not backbone_matches:
        return False, "Sphingosine backbone not found"
        
    # For each sphingosine backbone match found, verify that the nitrogen atom (the first atom in the pattern)
    # is acylated. This is done by checking that the nitrogen has a neighbor (other than the backbone carbon)
    # that is a carbonyl carbon (i.e. a carbon doubly bonded to an oxygen).
    acyl_found = False
    for match in backbone_matches:
        # The pattern is defined as: [0]=N, [1]=chiral C (attached to CO), [2]=chiral C (bearing OH)
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Get the indices of atoms connected to the nitrogen.
        neighbor_indices = [nbr.GetIdx() for nbr in n_atom.GetNeighbors()]
        # We know one neighbor (the backbone carbon) is part of the sphingosine; so check the other neighbor(s).
        for nbr_idx in neighbor_indices:
            # Skip the backbone neighbor if it is the one in our match (i.e. index match[1])
            if nbr_idx == match[1]:
                continue
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            # We expect the acyl group to be attached via an amide bond.
            # Thus the neighbor should be a carbon (atomic number 6).
            if nbr_atom.GetAtomicNum() == 6:
                # Check if this carbon (the acyl carbon) has an oxygen double-bonded
                for bond in nbr_atom.GetBonds():
                    # Look for an oxygen with a double bond to this carbon.
                    other = bond.GetOtherAtom(nbr_atom)
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        acyl_found = True
                        break
            if acyl_found:
                break
        if acyl_found:
            break

    if not acyl_found:
        return False, "Acyl group (via amide bond on the sphingosine nitrogen) not found"

    # If sphingosine backbone is present and its nitrogen is acylated, we classify as N-acylsphingosine.
    return True, "Molecule contains a sphingosine backbone with an N-linked acyl (fatty acid) group"