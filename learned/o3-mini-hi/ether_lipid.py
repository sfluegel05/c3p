"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether Lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol 
is bonded to an alkyl chain via an ether linkage (C-O-C) rather than the usual ester linkage.
This improved approach does not demand an exact match to a fixed glycerol SMARTS.
Instead, we require:
  (1) a molecular weight above a threshold (to avoid small molecules),
  (2) the presence of a C–O–C linkage that is not an ester (i.e. neither carbon is carbonyl),
  (3) one side of the ether linkage is “polar” – that is, it bears at least one additional oxygen substituent (hinting at a glycerol headgroup),
  (4) while the other side gives rise to a contiguous alkyl chain (of at least eight carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Rather than matching a fixed glycerol SMARTS, the function uses the following heuristic:
      - The molecule must have a molecular weight large enough for a lipid (here, >=300 Da).
      - At least one C–O–C (ether) linkage must be found in which neither carbon is carbonyl.
      - In that ether linkage one carbon must display polarity (i.e. have at least one oxygen neighbor other than
        the ether oxygen; a rough indicator for a glycerol-like headgroup) while the other carbon is attached to a long, 
        contiguous alkyl chain (with eight or more carbon atoms).
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Filter out very small molecules.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, "Molecular weight too low for a lipid"
        
    # Define a SMARTS pattern for an ether linkage: a carbon-oxygen-carbon unit.
    ether_smarts = "[#6]-O-[#6]"
    ether_query = Chem.MolFromSmarts(ether_smarts)
    ether_matches = mol.GetSubstructMatches(ether_query)
    if not ether_matches:
        return False, "No C-O-C ether linkage found"
        
    # Helper function: check whether a carbon atom is part of a carbonyl group.
    def is_carbonyl(carbon):
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:  # oxygen
                    return True
        return False
    
    # Helper: determine if an atom is an aliphatic (non-aromatic) carbon.
    def is_aliphatic_carbon(atom):
        return atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()
    
    # Helper: recursively determine the longest contiguous chain (by number of carbon atoms)
    # starting from a given carbon atom. We allow bonds that are single or double (to include unsaturated chains).
    def dfs_chain(atom, visited):
        length = 1
        visited.add(atom.GetIdx())
        for bond in atom.GetBonds():
            if bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and nbr.GetIdx() not in visited:
                    # Continue depth-first search from this neighbor.
                    new_length = 1 + dfs_chain(nbr, visited.copy())
                    if new_length > length:
                        length = new_length
        return length

    # For each ether linkage found, check our criteria.
    for match in ether_matches:
        # match is a tuple (carbon1_idx, oxygen_idx, carbon2_idx)
        c1 = mol.GetAtomWithIdx(match[0])
        o_atom = mol.GetAtomWithIdx(match[1])
        c2 = mol.GetAtomWithIdx(match[2])
        
        # Exclude ether bonds that are part of ester groups (one of the carbons is carbonyl).
        if is_carbonyl(c1) or is_carbonyl(c2):
            continue
        
        # We now consider the two sides of the ether linkage.
        # For each ordering, we will treat one carbon as "polar/glycerol-like" (must have at least one additional O neighbor)
        # and the other as the "alkyl chain" (which must link to a long alkyl chain, here defined as >=8 contiguous carbons).
        for polar_carbon, chain_carbon in [(c1, c2), (c2, c1)]:
            # Check if the polar candidate has at least one additional oxygen neighbor (other than the ether oxygen).
            has_extra_oxygen = False
            for nbr in polar_carbon.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != o_atom.GetIdx():
                    has_extra_oxygen = True
                    break
            if not has_extra_oxygen:
                continue  # This side does not appear polar enough
            
            # Now, for the candidate chain side, try to find a long contiguous alkyl chain.
            longest_chain = 0
            # Examine all neighbors of chain_carbon (except the oxygen already in the ether linkage).
            for nbr in chain_carbon.GetNeighbors():
                if nbr.GetIdx() == o_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    chain_len = dfs_chain(nbr, set())
                    if chain_len > longest_chain:
                        longest_chain = chain_len
            # If we have at least 8 connected carbon atoms then we consider this ether linkage as the one seen in an ether lipid.
            if longest_chain >= 8:
                return True, ("Molecule contains a glycerol-like moiety (carbon with extra oxygen substituent) "
                              "and an ether linkage to a long alkyl chain (chain length %d >=8)" % longest_chain)
    
    return False, "No suitable ether linkage attached to a glycerol-like headgroup and a long alkyl chain found"

# Example usage:
if __name__ == "__main__":
    # Provided test SMILES examples (only a few shown here for demonstration)
    smiles_examples = [
        "P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COCCCCCCCCCCCCCCCCCC)(O)(O)=O",
        "C[C@H]1CCC[C@@H](C)CCC[C@@H](C)CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCOC(CO)COCC[C@@H](C)CCC1",
        "CCCCCCCC\\C=C/CCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C",
        "C([C@@](CO/C=C\\CCCCCCCCCCCCCCCC)(O)[H])OP([O-])(=O)OCC[N+](C)(C)C"
    ]
    for s in smiles_examples:
        result, reason = is_ether_lipid(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")