"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N‐acylsphingosine (parent compounds of the ceramide family)

Definition: An N‐acylsphingosine is a molecule that contains a sphingosine backbone –
an acyclic chain featuring a secondary amine attached to two hydroxyl‐bearing carbons
(where one of these carbons is linked to a C=C bond) – and where the nitrogen is acylated
(via an amide bond to a fatty acyl group that itself has a long aliphatic chain).
"""

from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    Our approach is to look for:
      1. A sphingosine backbone defined by (ignoring stereochemistry):
         an N atom directly bonded to a carbon which is bonded to a CH2OH,
         then to a second carbon carrying an OH which in turn is connected to a C=C fragment.
      2. That the N atom is acylated; i.e. it has an additional neighbor (other than the backbone)
         that is a carbonyl carbon (bonded via a double bond to oxygen).
      3. That the acyl (fatty acid) group contains a long aliphatic chain (we require at least 6 connected carbons).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the N-acylsphingosine criteria, False otherwise.
        str: A message indicating the reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a sphingosine-backbone SMARTS.
    # This pattern (ignoring stereochemistry) requires:
    #   N - C(CO) - C(O) - C=C
    # which captures a N atom bound to a first carbon (with a CH2OH group CO),
    # then a second carbon with OH that is attached to an alkene.
    backbone_smarts = "N-C(CO)-C(O)C=C"
    backbone_pattern = Chem.MolFromSmarts(backbone_smarts)
    if backbone_pattern is None:
        return False, "Error in backbone SMARTS definition"
    
    backbone_matches = mol.GetSubstructMatches(backbone_pattern, useChirality=False)
    if not backbone_matches:
        return False, "Sphingosine backbone not found"
    
    # Helper function: perform DFS to determine longest chain of connected (aliphatic) carbons.
    def dfs_chain_length(atom, coming_from, visited):
        max_length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from:
                continue
            # Only count if the neighbor is carbon, non-aromatic.
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                if nbr.GetIdx() in visited:
                    continue
                visited.add(nbr.GetIdx())
                length = 1 + dfs_chain_length(nbr, atom.GetIdx(), visited)
                if length > max_length:
                    max_length = length
                visited.remove(nbr.GetIdx())
        return max_length

    # Set a minimum chain length for the fatty acyl group.
    min_chain_length = 6

    # Iterate over each sphingosine backbone match.
    for match in backbone_matches:
        # Check that the backbone atoms are acyclic (not in a ring)
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            continue  # skip if backbone is in a ring

        n_idx = match[0]  # the first atom is the nitrogen in our defined pattern.
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Now check that the N is acylated.
        acyl_found = False
        for nbr in n_atom.GetNeighbors():
            # We expect one neighbor to be part of the backbone.
            if nbr.GetIdx() in match:
                continue
            # Look for a carbon neighbor that is a potential carbonyl carbon.
            if nbr.GetAtomicNum() == 6:
                # Check for a double bond to oxygen from this neighbor.
                carbonyl_flag = False
                for bond in nbr.GetBonds():
                    # Look for double bond to oxygen.
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            carbonyl_flag = True
                            break
                if not carbonyl_flag:
                    continue

                # We have now a candidate acyl carbon. Check that it carries a long aliphatic chain.
                # To do so, we search from this carbon (excluding the bond back to the N)
                visited = set([nbr.GetIdx()])
                chain_length = dfs_chain_length(nbr, n_idx, visited)
                # Include the carbonyl carbon itself as part of the chain.
                total_chain = 1 + chain_length
                if total_chain >= min_chain_length:
                    acyl_found = True
                    break

        if not acyl_found:
            continue  # try next backbone match
        
        # If we reached here then we have a sphingosine backbone, the N is acylated with a valid carbonyl,
        # and the acyl chain is long enough.
        return True, "Molecule contains a sphingosine backbone with an N-linked acyl (fatty acid) group and a long aliphatic chain"
    
    return False, "No valid N-acylsphingosine backbone with proper acylation and long aliphatic chain found"


# Example usage:
if __name__ == "__main__":
    # Test one of the provided examples: N-2-hydroxylignoceroylsphingosine
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_N_acylsphingosine(test_smiles)
    print(result, reason)