"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
Definition: Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.
In our algorithm, we attempt to locate a glycerol backbone (a contiguous three-carbon fragment)
in which each carbon bears one oxygen substituent. In a free glycerol all three oxygens would be –OH;
here we expect that exactly one oxygen is substituted (no attached H) while the remaining two are free –OH.
Additionally, the substituent group (R) is expected to be a carbon chain of a minimal length.
If these conditions are met, we return True with a reason; otherwise, False.
"""

from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    It looks for a glycerol backbone (three contiguous sp3 carbons that in a “free” glycerol 
    would each be attached to one oxygen) where exactly one oxygen substituent is not a hydroxyl and 
    the other two are free hydroxyls. In addition, the substituent group is required to have a minimum
    carbon chain length (here set to 4) to rule out very short substitutions.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoradylglycerol, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function to compute the longest connected carbon chain starting from a given carbon atom.
    def longest_carbon_chain(atom, parent_idx, visited):
        # Only count carbon atoms
        if atom.GetAtomicNum() != 6:
            return 0
        visited.add(atom.GetIdx())
        max_length = 1  # count current atom
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == parent_idx: 
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + longest_carbon_chain(nbr, atom.GetIdx(), visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    # We now attempt to find a glycerol-like backbone,
    # i.e. a linear chain of three connected carbons (C1 - C2 - C3)
    # Expectation: Each carbon has one oxygen substituent that is not part of the chain.
    # In free glycerol, all three oxygens would be -OH. For a monoradylglycerol exactly one oxygen shall be "substituted"
    # (no attached H) while the other two are free hydroxyls (with at least one hydrogen).
    
    # Loop over atoms to find candidate middle carbons (C2) that have at least two sp3 carbon neighbors.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Only consider sp3 carbons (for a typical glycerol backbone)
        if atom.GetHybridization().name != "SP3":
            continue
        # Collect neighboring carbons that are sp3
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization().name=="SP3"]
        if len(carbon_neighbors) < 2:
            continue
        # For each distinct pair of neighboring carbons, try a candidate chain: (nbr1, current, nbr2).
        for i in range(len(carbon_neighbors)):
            for j in range(i+1, len(carbon_neighbors)):
                c1 = carbon_neighbors[i]
                c3 = carbon_neighbors[j]
                # Avoid cases where the two neighbors are directly connected (would be a ring) 
                if c1.IsInRing() and c3.IsInRing():
                    # skip if they are connected by a bond (not expected for a linear glycerol backbone)
                    if c1.HasSubstructMatch(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmarts(f"[#{c1.GetAtomicNum()}]~[#{c3.GetAtomicNum()}]")))):
                        pass  # not a good check; instead we simply check bond existence:
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()):
                    continue
                    
                # Our candidate backbone is the three carbons: c1, atom (as c2), c3.
                backbone = [c1, atom, c3]
                # For each backbone carbon, find oxygen neighbors that are not part of the backbone.
                oxygens = []
                valid_backbone = True
                for carbon in backbone:
                    # Get oxygen neighbors that are not in the candidate chain
                    o_neighbors = [nbr for nbr in carbon.GetNeighbors() 
                                   if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [a.GetIdx() for a in backbone]]
                    # We expect one oxygen substituent per carbon in glycerol (if fully substituted as in glycerol triol)
                    if len(o_neighbors) != 1:
                        valid_backbone = False
                        break
                    oxygens.append(o_neighbors[0])
                    
                if not valid_backbone:
                    continue
                    
                # Now check whether exactly two of the three oxygens are free hydroxyl groups
                # and one is substituted (no attached hydrogen).
                free_OH_count = 0
                substituted_idx = -1
                for idx, o_atom in enumerate(oxygens):
                    # use GetTotalNumHs() to check for H (includes implicit)
                    numHs = o_atom.GetTotalNumHs()
                    if numHs > 0:
                        free_OH_count += 1
                    else:
                        substituted_idx = idx
                if free_OH_count != 2 or substituted_idx < 0:
                    continue

                # Next, verify that the substituent group (attached via the oxygen that is not free)
                # is a carbon-based chain of a minimum length threshold.
                # Get the substituted oxygen and then its neighbor that is not the backbone carbon.
                sub_oxygen = oxygens[substituted_idx]
                # Its neighbors: one is the glycerol carbon
                sub_neighbors = [nbr for nbr in sub_oxygen.GetNeighbors() if nbr.GetIdx() not in [a.GetIdx() for a in backbone]]
                if not sub_neighbors:
                    continue
                # We expect the substituent to be carbon-based.
                sub_atom = sub_neighbors[0]
                if sub_atom.GetAtomicNum() != 6:
                    continue
                chain_length = longest_carbon_chain(sub_atom, sub_oxygen.GetIdx(), set())
                # Here we choose a minimum chain length of 4 carbons as a heuristic.
                if chain_length < 4:
                    continue

                # If we reach here, we found a candidate glycerol backbone with exactly one substitution.
                reason = (f"Found glycerol backbone with substitution at position {substituted_idx+1}; "
                          f"chain length of substituent = {chain_length}, free hydroxyls = 2.")
                return True, reason

    # If no candidate backbone is found, reject.
    return False, "No glycerol backbone with a single substituted position (and two free OH groups) was found."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "O=C(OC[C@@H](O)CO)CCCCCCC/C=C(\\CCCCCCCC)/C",  # 2,3-dihydroxypropyl (Z)-10-methyloctadec-9-enoate
        "CCCCCCCC(=O)OCC(O)CO",                         # 1-monooctanoylglycerol
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO",           # 3-stearoyl-sn-glycerol
        "O(C(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO",  # MG(22:2(13Z,16Z)/0:0/0:0)
        "O(CC(O)CO)C(=O)CC",                            # Glycerol 1-propanoate
    ]
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\n Result: {result} | Reason: {reason}\n")