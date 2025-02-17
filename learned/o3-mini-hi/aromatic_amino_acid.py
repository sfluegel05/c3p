"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
Definition: An amino acid whose structure includes an aromatic ring.
The molecule must contain an integrated (free) alpha–amino acid backbone,
i.e. an sp3 alpha–carbon with a free (non‐amidated) amino group (NH2 or NH3+)
and a free carboxyl group (COOH, COO–), and the side chain attached to the alpha–carbon
must contain an aromatic ring.
"""

from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines whether a molecule is an aromatic amino acid.
    A valid aromatic amino acid must have a free alpha–amino acid backbone
    (an sp3 alpha–carbon attached to a free (non‐amidic) amino group and a free –COOH group)
    and the side chain must include an aromatic ring.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies, False otherwise.
        str: An explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens to ensure proper matching.
    mol = Chem.AddHs(mol)
    
    # Iterate over potential alpha–carbon atoms. We require the alpha–carbon to be sp3,
    # non‐aromatic and to have exactly three neighbors (one for the amine, one for the carboxyl,
    # one for the side chain). (Note: Glycine will be skipped since its side chain is just H.)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetIsAromatic():
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 3:
            continue  # not the typical connectivity for an isolated amino acid
        # We will try to identify three groups:
        #   1. A free (non‐amidated) amino group.
        #   2. A carboxyl group (COOH or COO–) that is not part of a peptide bond.
        #   3. A side chain.
        amine_atom = None
        carboxyl_atom = None
        side_chain_atom = None
        
        for nbr in neighbors:
            # Check for non‐amidated amino group: a nitrogen with at least two attached hydrogens.
            if nbr.GetAtomicNum() == 7:
                # Count explicit hydrogens on the nitrogen.
                h_count = sum(1 for a in nbr.GetNeighbors() if a.GetAtomicNum() == 1)
                # In a free amino group we expect at least 2 H's.
                if h_count >= 2:
                    amine_atom = nbr
                    continue
            # Check for carboxyl group: a carbon (should be sp2) that has a double bond to an O 
            # and a single-bonded O with at least one hydrogen, and is not further amidated.
            if nbr.GetAtomicNum() == 6:
                # Preliminary: carboxyl carbon should have at least 3 neighbors (one being the alpha carbon).
                if len(nbr.GetNeighbors()) < 2:
                    continue
                oxy_double = False
                oxy_hydroxyl = False
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        # Check for carbonyl oxygen.
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            oxy_double = True
                        # Check for hydroxyl oxygen (SINGLE bond and should have a hydrogen).
                        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            if any(neigh.GetAtomicNum() == 1 for neigh in subnbr.GetNeighbors()):
                                oxy_hydroxyl = True
                if oxy_double and oxy_hydroxyl:
                    carboxyl_atom = nbr
                    continue
            # Anything else will be treated as the candidate side chain.
            side_chain_atom = nbr
        
        # Check if we found all three parts.
        if not (amine_atom and carboxyl_atom and side_chain_atom):
            continue  # this alpha–carbon candidate is missing one of the groups
        
        # At this point we have a candidate alpha–carbon with:
        #  - one neighbor as a free amine,
        #  - one neighbor as a carboxyl group,
        #  - one neighbor as the side chain.
        # Now we ensure that the side chain (when taken as the fragment, avoiding the backbone atoms)
        # contains at least one aromatic atom.
        
        # We perform a depth-first search starting from the side_chain_atom.
        # Do not traverse back to any of the backbone atoms: the alpha–carbon, the amine, or carboxyl.
        backbone_indices = {atom.GetIdx(), amine_atom.GetIdx(), carboxyl_atom.GetIdx()}
        visited = set()
        stack = [side_chain_atom.GetIdx()]
        aromatic_found = False
        while stack:
            current_idx = stack.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            if current_atom.GetIsAromatic():
                aromatic_found = True
                break
            # Traverse neighbors if they are not in the backbone.
            for nbr in current_atom.GetNeighbors():
                if nbr.GetIdx() not in backbone_indices and nbr.GetIdx() not in visited:
                    stack.append(nbr.GetIdx())
        if aromatic_found:
            return True, "Contains integrated alpha–amino acid backbone with an aromatic side chain"
        else:
            return False, "Alpha–amino acid backbone found but no aromatic ring connected in the side chain"
    
    # If no candidate alpha–carbon satisfying conditions is found, then we did not match the backbone.
    return False, "No integrated free alpha–amino acid backbone found"

# Example usage:
if __name__ == "__main__":
    # Test with D-histidine (one of the true positives)
    test_smiles = "N[C@H](Cc1c[nH]cn1)C(O)=O"
    result, reason = is_aromatic_amino_acid(test_smiles)
    print(result, reason)