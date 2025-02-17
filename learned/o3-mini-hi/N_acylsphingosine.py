"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N‐acylsphingosine (parent compounds of the ceramide family)
Definition: Molecules that contain a sphingosine backbone – that is, an acyclic chain 
featuring a secondary amine bonded to two hydroxyl‐bearing carbons – where the nitrogen 
is acylated (i.e. forms an amide bond to a fatty acyl group). In addition, a bona fide 
N‐acylsphingosine should have a long aliphatic (lipid) chain.
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is defined as a sphingosine backbone (an acyclic chain featuring
    a secondary amine attached to two hydroxylated carbons) in which the nitrogen is 
    acylated via an amide bond to a fatty acyl group. In addition, a long aliphatic 
    chain (a proxy for lipid character) is expected.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule matches the N-acylsphingosine criteria, False otherwise.
        str: A message indicating the reason for the decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Use a relaxed backbone SMARTS pattern for sphingosine.
    # We want to capture a nitrogen attached to two carbons bearing hydroxyl groups.
    # We ignore chirality details by not encoding stereochemistry.
    # This pattern is: N-C(CO)(O)
    backbone_pattern = Chem.MolFromSmarts("N-C(CO)(O)")
    if backbone_pattern is None:
        return False, "Error in backbone SMARTS definition"
        
    # Search for matches of the sphingosine backbone.
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if not backbone_matches:
        return False, "Sphingosine backbone not found"
        
    # For each matching backbone, check extra criteria.
    for match in backbone_matches:
        # 1. Ensure that the three backbone atoms (N and the two carbons) are acyclic.
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            continue  # Skip matches that occur in rings (likely sugars or cyclic fragments)

        # 2. Check that the nitrogen (first atom in the pattern) is acylated.
        #    In other words, aside from the backbone connection, the nitrogen should be 
        #    bonded to a carbon that is part of an amide group (i.e. that carbon is double-bonded to oxygen).
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        acyl_found = False
        for nbr in n_atom.GetNeighbors():
            # Skip the backbone neighbor (it should be one of the atoms in the backbone match)
            if nbr.GetIdx() in match:
                continue
            # Look for a carbon neighbor
            if nbr.GetAtomicNum() == 6:
                # Check if any bond from that carbon is a double bond to oxygen (carbonyl).
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            # To avoid matching groups that are embedded in rings,
                            # require the carbonyl carbon is not in a ring.
                            if not nbr.IsInRing():
                                acyl_found = True
                                break
                if acyl_found:
                    break

        if not acyl_found:
            continue  # Try the next backbone match
        
        # 3. Check that the molecule contains a long aliphatic chain.
        # This is a proxy for a fatty acid or long hydrocarbon tail.
        # Here we look for a substructure of 6 consecutive carbon atoms.
        long_chain_pattern = Chem.MolFromSmarts("CCCCCC")
        if long_chain_pattern is None:
            return False, "Error defining long chain SMARTS"
        if not mol.HasSubstructMatch(long_chain_pattern):
            continue  # If there is no sufficiently long chain, skip this match
        
        # If all criteria (acyclic backbone, acylation, long aliphatic chain) are met,
        # then classify as an N-acylsphingosine.
        return True, "Molecule contains an acyclic sphingosine backbone with an N-linked acyl (fatty acid) group and a long aliphatic chain"
    
    # If no match fulfills all criteria, the molecule is not classified as N-acylsphingosine.
    return False, "No valid N-acylsphingosine backbone with proper acylation and long aliphatic chain found"

# Example usage (these can be replaced or removed in production):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC"  # one example from the list
    result, reason = is_N_acylsphingosine(test_smiles)
    print(result, reason)