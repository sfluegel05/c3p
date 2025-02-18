"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (free or as an acyl chain) that contains at least one C=C double bond.
A fatty acid is defined as a linear, unbranched chain (attached via a free acid or acyl ester “handle”)
of at least 6 contiguous carbon atoms that has at least one carbon–carbon double bond.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

# Helper function to check if an atom is a carbon and is not in a ring.
def is_linear_carbon(atom):
    return atom.GetAtomicNum() == 6 and not atom.IsInRing()

# Helper function: starting from a carbon atom that is presumed to be terminal (has one neighbor aside from the handle),
# walk along the chain (always following the unique next carbon) and return:
#   length: number of carbons in this linear chain
#   has_double: whether at least one bond in the chain is a C=C double bond.
def get_linear_chain_info(atom, mol, coming_from_idx):
    # We count the current atom
    length = 1
    has_double = False
    current_atom = atom
    current_idx = atom.GetIdx()
    while True:
        # Get neighbors that are carbon and not the atom we came from.
        nbrs = [nbr for nbr in current_atom.GetNeighbors() 
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != coming_from_idx and not nbr.IsInRing()]
        # In a truly linear (unbranched) fatty acyl chain the current carbon should have at most one next carbon.
        if len(nbrs) != 1:
            break
        next_atom = nbrs[0]
        bond = mol.GetBondBetweenAtoms(current_idx, next_atom.GetIdx())
        if bond is None:
            break
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            has_double = True
        # Update for next step
        coming_from_idx = current_idx
        current_atom = next_atom
        current_idx = current_atom.GetIdx()
        length += 1
    return length, has_double

# Main function: classify if the molecule (given as a SMILES string) is an olefinic fatty acid.
def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    
    The rules:
      1. Look for either a free carboxylic acid group (defined by pattern "C(=O)[O;H1,O-]")
         or an acyl ester motif ("C(=O)O[C]").
      2. For each match, identify the “handle” carbon (the one attached to the C=O group) that is expected
         to be the start of a fatty acyl chain.
      3. The handle must be terminal (i.e. it has no extra carbon neighbors except the one connecting 
         to the acid/ester carbonyl).
      4. From this handle, “walk” along a strictly linear carbon chain. Require that the chain includes 
         at least 6 carbons and that at least one bond in the chain is a C=C double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an olefinic fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    MIN_CHAIN_LENGTH = 6  # Minimum number of carbon atoms in the acyl chain

    # Quick global screening: there must be a C=C somewhere in the molecule.
    dbl_bond_smarts = Chem.MolFromSmarts("[C]=[C]")
    if not mol.HasSubstructMatch(dbl_bond_smarts):
        return False, "No carbon–carbon double bond (C=C) found in molecule"
    
    reasons = []
    # First, try free acid group: expect a free acid group to be defined by "C(=O)[O;H1,O-]"
    free_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    free_matches = mol.GetSubstructMatches(free_acid_pattern)
    if free_matches:
        for match in free_matches:
            # In the free acid SMARTS, match[0] is the carbonyl carbon.
            acid_carbon = mol.GetAtomWithIdx(match[0])
            # For a free fatty acid, the acid carbon should have exactly one neighboring carbon (the alpha carbon).
            carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                reasons.append("Free acid group is not terminal (acid carbon has >1 carbon neighbor)")
                continue
            alpha_c = carbon_neighbors[0]
            # Ensure the handle is a linear (non-branched) carbon.
            # The alpha carbon should have only one carbon neighbor other than the acid carbon.
            alpha_neigh = [nbr for nbr in alpha_c.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != acid_carbon.GetIdx()]
            if len(alpha_neigh) != 1:
                reasons.append("Alpha carbon for free acid is branched")
                continue
            # Walk down the chain starting from the alpha carbon.
            chain_length, has_double = get_linear_chain_info(alpha_c, mol, acid_carbon.GetIdx())
            if chain_length < MIN_CHAIN_LENGTH:
                reasons.append(f"Free acid chain too short (length {chain_length} < {MIN_CHAIN_LENGTH})")
                continue
            if not has_double:
                reasons.append("Free acid chain does not contain a C=C double bond")
                continue
            return True, "Contains a fatty acyl chain (via free acid) with sufficient length and a C=C double bond."
    
    # Next, try acyl ester motif: defined as the three-atom substring "C(=O)O[C]"
    acyl_ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if ester_matches:
        for match in ester_matches:
            # In this SMARTS, match[0] is the carbonyl carbon, match[1] is the ester oxygen and match[2] is the handle.
            if len(match) < 3:
                continue
            handle_atom = mol.GetAtomWithIdx(match[2])
            # For acyl esters, we similarly require that the handle is terminal.
            handle_neigh = [nbr for nbr in handle_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != match[1]]
            if len(handle_neigh) != 1:
                reasons.append("Acyl ester handle is branched")
                continue
            chain_length, has_double = get_linear_chain_info(handle_atom, mol, match[1])
            if chain_length < MIN_CHAIN_LENGTH:
                reasons.append(f"Acyl ester chain too short (length {chain_length} < {MIN_CHAIN_LENGTH})")
                continue
            if not has_double:
                reasons.append("Acyl ester chain does not contain a C=C double bond")
                continue
            return True, "Contains a fatty acyl chain (via acyl ester) with sufficient length and a C=C double bond."
    
    # If nothing qualified, report the reasons encountered (if any) or a default message.
    if reasons:
        return False, " ; ".join(reasons)
    else:
        return False, "No free acid or acyl ester substructure with a qualifying fatty acyl chain found."

# Example usage (for testing; remove or comment out during integration).
if __name__ == "__main__":
    test_smiles_list = [
        "O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O",  # Avenoleic acid (should be True)
        "O=C(CCC/C=C\\C/C=C\\CC(/C=C/C(C(CCCCC)O)O)O)O",  # 11,14,15-trihydroxy-(5Z,8Z,12E)-icosatrienoic acid (True)
        "CC\\C=C/C\\C=C/CCC\\C=C\\C=C\\C=C/CCCC(O)=O",  # (5Z,7E,9E,14Z,17Z)-icosapentaenoic acid (True)
        "O(C(C[N+](C)(C)C)CC([O-])=O)C(=O)CC/C=C/C/C=C/C/C=C/C\\C=C\\CC(O)/C=C/C=C/CC",  # False positive example from phospholipid-like structure
        "O=C1N(C(=O)C=2N=CN=NC2N1C)",  # Fervenulin (should be False: no C=C present)
    ]
    
    for s in test_smiles_list:
        res, msg = is_olefinic_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*60}")