"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: An aliphatic alcohol, defined as 'An alcohol derived from an aliphatic compound.'
Improvement rationale:
 - First disqualify molecules that contain free carboxylic acid/carboxylate groups.
 - For each –OH group, require that it is attached to a saturated (sp³), non‐aromatic carbon.
 - Then, from that candidate carbon (ignoring the –OH oxygen) perform a DFS search on contiguous 
   chain carbons that are acyclic and non‐aromatic. This DFS now allows a “tolerance” of up to 
   two isolated sp² centers (if detected) but disallows carbons that are carbonyl-like.
 - If any adjoining chain has at least 6 carbons (not counting the candidate carbon with -OH),
   then the molecule is classified as an aliphatic alcohol.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    Here, an aliphatic alcohol is defined as having at least one hydroxyl (-OH) group attached 
    to a saturated (sp³), non‐aromatic carbon. Furthermore, that candidate carbon must be connected 
    to a contiguous, open-chain aliphatic (acyclic) region (which may include one or two isolated sp2 
    centers) of at least 6 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Disqualify molecules with free carboxylic acid/carboxylate functionality.
    acid_smarts = [Chem.MolFromSmarts("[CX3](=O)[O;H1]"), 
                   Chem.MolFromSmarts("[CX3](=O)[O-]")]
    for pattern in acid_smarts:
        if mol.HasSubstructMatch(pattern):
            return False, "Molecule contains carboxylic acid/carboxylate functionality"
    
    # --- Helper function to check if a carbon atom can be part of the chain.
    # Returns:
    #   "sp3" if the atom is sp3 (always acceptable),
    #   "sp2" if the atom is sp2 (allowed but only isolated),
    #   False if not allowed.
    def allowed_chain_atom(atom):
        if atom.GetAtomicNum() != 6:
            return False
        # Must not be aromatic or in a ring (acyclic chain).
        if atom.GetIsAromatic() or atom.IsInRing():
            return False
        # Exclude carbons that are carbonyl-like (i.e. have a double bond to oxygen).
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return False
        hyb = atom.GetHybridization()
        if hyb == Chem.rdchem.HybridizationType.SP3:
            return "sp3"
        elif hyb == Chem.rdchem.HybridizationType.SP2:
            return "sp2"
        return False

    # --- DFS search: starting at a chain atom, count the maximum contiguous chain length.
    # We allow a penalty counter for sp2 atoms (max allowed is 2 per path).
    def dfs_chain(atom, visited, sp2_count):
        visited.add(atom.GetIdx())
        max_length = 1  # count this atom
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6 or nbr.GetIdx() in visited:
                continue
            typ = allowed_chain_atom(nbr)
            if not typ:
                continue
            new_sp2_count = sp2_count
            if typ == "sp2":
                new_sp2_count += 1
                if new_sp2_count > 2:  # exceed tolerance
                    continue
            branch_length = dfs_chain(nbr, visited.copy(), new_sp2_count)
            if branch_length + 1 > max_length:
                max_length = branch_length + 1
        return max_length

    # --- Iterate over oxygen atoms that may belong to a hydroxyl (-OH) group.
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum() != 8:
            continue
        # Check that oxygen is a free hydroxyl (has at least one hydrogen present).
        if oxygen.GetTotalNumHs() < 1:
            continue
        # Check each neighbor of oxygen for a candidate attachment carbon.
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # Candidate carbon must be sp3 and non-aromatic.
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or nbr.GetIsAromatic():
                continue
            # From this candidate carbon, search its other neighbors (excluding the hydroxyl oxygen)
            # to find an adjoining acyclic chain.
            for subnbr in nbr.GetNeighbors():
                if subnbr.GetIdx() == oxygen.GetIdx():
                    continue  # do not go back to the -OH oxygen
                typ = allowed_chain_atom(subnbr)
                if not typ:
                    continue
                # Perform a DFS starting from subnbr to determine contiguous chain length.
                chain_length = dfs_chain(subnbr, set(), 0)
                if chain_length >= 6:
                    return True, ("Found -OH group attached to a sp3, non-aromatic carbon that is connected "
                                  "to an acyclic aliphatic chain (chain length = {}).".format(chain_length))
    return False, ("No qualifying -OH group found attached to a sp3, non-aromatic carbon with an adjoining "
                   "acyclic aliphatic chain of at least 6 carbons.")

# (Optional) Testing examples.
if __name__ == "__main__":
    test_smiles = [
        # Examples expected to classify as aliphatic alcohols (true positives):
        "CCCCCCC(C)O",  # octan-2-ol (chain length = 6)
        "O=C(OC)/C=C/CC(O)CCCCCC(O)C",  # Cladosporester A (chain length = 9)
        "CCCCCC(O)CC(=O)CCc1ccc(O)c(O)c1",  # 1-(3,4-Dihydroxyphenyl)-5-hydroxy-3-decanone (chain length = 7)
        "OCCCCCC/C=C/CCCC",  # 7E-dodecen-1-ol (chain length = 11)
        "OCCCCCCCCCCC1=CC=CC=C1",  # 10-Phenyl-1-decanol (chain length = 9)
        "CCCCCCCCCCCCCCCCCCCCCCO",  # tricosan-1-ol (chain length = 22)
        # False positive candidate: contains carboxylate functionality.
        "CCCCCCCCCCCCCCCCCC(O)C([O-])=O",  # 2-hydroxyarachidate
    ]
    
    for sm in test_smiles:
        classification, reason = is_aliphatic_alcohol(sm)
        print("SMILES:", sm)
        print("Classification:", classification)
        print("Reason:", reason)
        print("-" * 60)