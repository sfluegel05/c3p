"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (free acid or as an acyl chain) that contains at least one C=C double bond.
A fatty acid here is defined as a linear (i.e. strictly unbranched) chain of at least 6 contiguous carbons
attached either to a free acid group or as an acyl ester chain; and must contain at least one carbon–carbon double bond.
Improvements over the previous method:
  • Instead of a DFS that aborts as soon as >1 candidate neighbor is encountered, we now recursively
    extract the chain length and check for branch points. If any branch is found the candidate is disqualified.
  • We “start” the chain from a candidate handle atom (alpha carbon from the acid, or acyl ester handle)
    and then follow the unique contiguous carbon path (ignoring ring atoms) while checking for at least one double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

MIN_CHAIN_LENGTH = 6  # at least 6 contiguous carbons

def get_linear_chain(mol, current_idx, prev_idx):
    """
    Recursively determine the length of a linear (unbranched) chain starting from current_idx.
    Also checks if at least one traversed bond is a double bond (C=C).
    For our strict fatty acid definition, if a branch is encountered (more than one candidate neighbor),
    the chain is disqualified.
    
    Args:
       mol: RDKit molecule.
       current_idx: current atom index (should be carbon).
       prev_idx: the previous atom's index from which we came.
       
    Returns:
       (chain_length, found_double, branched) where:
         chain_length (int): number of carbon atoms in the linear chain (including the starting atom).
         found_double (bool): True if at least one double bond was encountered along the chain.
         branched (bool): True if a branch (more than one candidate) was found.
    """
    current_atom = mol.GetAtomWithIdx(current_idx)
    # Consider only neighbor carbons that are not the previous and are not in rings.
    candidates = []
    for nbr in current_atom.GetNeighbors():
        if nbr.GetIdx() == prev_idx:
            continue
        if nbr.GetAtomicNum() != 6:
            continue
        if nbr.IsInRing():
            continue
        candidates.append(nbr)
    
    if len(candidates) == 0:
        # Terminal carbon: chain of length 1 with no double bonds encountered in this step.
        return (1, False, False)
    if len(candidates) > 1:
        # Branching has occurred; we disqualify this candidate chain.
        return (0, False, True)
    
    # Exactly one candidate: follow that bond.
    next_atom = candidates[0]
    bond = mol.GetBondBetweenAtoms(current_idx, next_atom.GetIdx())
    bond_is_double = (bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE)
    # Recursively get the chain info from the next atom.
    child_length, child_found_double, child_branched = get_linear_chain(mol, next_atom.GetIdx(), current_idx)
    if child_branched:
        return (0, False, True)
    total_length = 1 + child_length
    found_double = bond_is_double or child_found_double
    return (total_length, found_double, False)

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as an olefinic fatty acid.
    
    Rules:
       1. There must be at least one C=C double bond in the molecule.
       2. Look for a free acid motif (C(=O)[O;H1,O-]) or an acyl ester motif (C(=O)O[C]).
       3. For a free acid, take as the candidate handle one of the alpha carbon neighbors (non‐ring carbon).
          For an acyl ester, use the handle as the carbon attached to the ester oxygen.
       4. From the candidate handle, examine the contiguous linear (unbranched) carbon chain
          (using get_linear_chain) that must have at least MIN_CHAIN_LENGTH carbons and at least one C=C double bond.
    
    Returns:
         (bool, str): (True, reason) if the molecule qualifies;
                      (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Quick global check: there must be at least one carbon–carbon double bond anywhere.
    dbl_bond_smarts = Chem.MolFromSmarts("[#6]=[#6]")
    if not mol.HasSubstructMatch(dbl_bond_smarts):
        return False, "No carbon–carbon double bond (C=C) found in molecule"
    
    reasons = []
    
    # FIRST: search for free acid motif.
    free_acid_smarts = "C(=O)[O;H1,O-]"  # matches free carboxylic acid
    free_acid_pattern = Chem.MolFromSmarts(free_acid_smarts)
    free_matches = mol.GetSubstructMatches(free_acid_pattern)
    
    if free_matches:
        for match in free_matches:
            # In our free acid motif, match[0] is the carbonyl carbon.
            acid_carbon = mol.GetAtomWithIdx(match[0])
            # Look for alpha carbon neighbor(s): carbons (not in ring) bonded to the acid carbon.
            alpha_candidates = [nbr for nbr in acid_carbon.GetNeighbors() 
                                if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
            if not alpha_candidates:
                reasons.append("Free acid group with no eligible alpha carbon found")
                continue
            for alpha in alpha_candidates:
                # Walk from the alpha carbon along the chain.
                chain_length, found_double, branched = get_linear_chain(mol, alpha.GetIdx(), acid_carbon.GetIdx())
                if branched:
                    # Disqualify chain if branching is present.
                    reasons.append("Free acid candidate chain is branched")
                    continue
                if chain_length < MIN_CHAIN_LENGTH:
                    reasons.append(f"Free acid candidate chain too short (length {chain_length})")
                    continue
                if not found_double:
                    reasons.append("Free acid candidate chain lacks a C=C double bond")
                    continue
                return True, ("Contains a fatty acyl chain (via free acid) with sufficient linear length "
                             "and a C=C double bond.")
        reasons.append("None of the free acid chains qualified (either too short, branched, or lacking a C=C double bond)")
    else:
        reasons.append("No free acid substructure (C(=O)[O;H1,O-]) found")
    
    # SECOND: search for acyl ester motifs.
    acyl_ester_smarts = "C(=O)O[C]"  # pattern for an acyl ester
    acyl_ester_pattern = Chem.MolFromSmarts(acyl_ester_smarts)
    ester_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    
    if ester_matches:
        for match in ester_matches:
            if len(match) < 3:
                continue
            # In the matched pattern: match[0] = carbonyl carbon, match[1] = ester oxygen, match[2] = handle carbon.
            handle_atom = mol.GetAtomWithIdx(match[2])
            # Walk from the handle along the chain.
            chain_length, found_double, branched = get_linear_chain(mol, handle_atom.GetIdx(), match[1])
            if branched:
                reasons.append("Acyl ester candidate chain is branched")
                continue
            if chain_length < MIN_CHAIN_LENGTH:
                reasons.append(f"Acyl ester candidate chain too short (length {chain_length})")
                continue
            if not found_double:
                reasons.append("Acyl ester candidate chain lacks a C=C double bond")
                continue
            return True, ("Contains a fatty acyl chain (via acyl ester) with sufficient linear length "
                         "and a C=C double bond.")
        reasons.append("Acyl ester chains did not qualify (either too short, branched, or lacking a C=C double bond)")
    else:
        reasons.append("No acyl ester substructure (C(=O)O[C]) found")
    
    return False, " ; ".join(reasons)

# Example usage and testing:
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O",  # Avenoleic acid
        "O=C(CCC/C=C\\C/C=C\\CC(/C=C/C(C(CCCCC)O)O)O)O",  # 11,14,15-trihydroxy-(5Z,8Z,12E)-icosatrienoic acid
        "CC\\C=C/C\\C=C/CCC\\C=C\\C=C\\C=C/CCCC(O)=O",  # (5Z,7E,9E,14Z,17Z)-icosapentaenoic acid
        "[H]C(\\C=C/CCCCCCCC(O)=O)=C1O[C@H]1C\\C=C/CC",  # (9Z,13S,15Z)-12,13-epoxyoctadeca-9,11,15-trienoic acid
        "OC(CCCCC)\\C=C\\C=C\\CCCCCCCC(O)=O",  # alpha-Artemisic acid
        "CCCCCCCCCC\\C=C/CCCCCCCC(O)=O",  # gadoleic acid (simulated longer chain)
        # False positives (should be classified as False):
        "O=C(O)/C=C/C=C/[C@@H]([C@@H]1O[C@]2(O[C@@H]([C@@H](CC)[C@@H](C2)O)C)[C@@H](C)[C@@H]([C@@H]1C)O)C",  
          # Pteridic acid F – has a free acid but its acyl candidate chain is not strictly linear
        # False negatives (should be classified as True):
        "OC[C@H](COP(OCC[N+](C)(C)C)(=O)[O-])OC(=O)CCCCCCC/C=C\\CCCCCC",  
          # 2-[(9Z)-hexadecenoyl]-sn-glycero-3-phosphocholine – acyl ester handle chain hopefully qualifies now.
    ]
    
    for sm in test_smiles:
        result, reason = is_olefinic_fatty_acid(sm)
        print("SMILES:", sm)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)