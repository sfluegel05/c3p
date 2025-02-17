"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (free or as an acyl chain) that contains at least one C=C double bond.
A fatty acid is identified by the presence of either a carboxylic acid group (free acid) or an acyl ester group
and by having a long, predominantly aliphatic carbon chain. For this classifier we require the chain length
(from the acid/ester “handle”) to be at least 6 carbon atoms and to contain at least one C=C bond.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

# Helper function: Find the longest contiguous chain (list of atom indices) that consists only of carbon atoms.
# We use DFS starting from a given atom and return the longest carbon path.
def longest_carbon_chain(atom, mol, visited=None):
    if visited is None:
        visited = set()
    curr_idx = atom.GetIdx()
    visited.add(curr_idx)
    best_path = [curr_idx]
    # Look at neighbors that are carbons and not yet visited.
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() != 6:
            continue
        nbr_idx = nbr.GetIdx()
        if nbr_idx in visited:
            continue
        # Recurse with a copy of visited.
        path = [curr_idx] + longest_carbon_chain(nbr, mol, visited.copy())
        if len(path) > len(best_path):
            best_path = path
    return best_path

# Check if there is at least one double bond (C=C) in the given chain.
def chain_has_double_bond(chain_atom_indices, mol):
    # iterate over consecutive pairs in the chain (in order)
    # Note: the DFS path may not be “linear” if branched, but our DFS returns one linear path
    for i in range(len(chain_atom_indices) - 1):
        a_idx = chain_atom_indices[i]
        b_idx = chain_atom_indices[i + 1]
        bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
        if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
            return True
    return False

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    
    It now does more refined checks:
    1. Look for either a free carboxylic acid group (C(=O)[O;H1,O-]) or an acyl ester motif (C(=O)O[C]).
    2. For each match, take the “handle” carbon (alpha carbon) connected to the acid (or ester) group.
    3. From the starting carbon, compute the longest contiguous (carbon-only) chain.
    4. Require that this chain is at least 6 carbons long.
    5. Require that at least one bond in the chain is a C=C double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an olefinic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Quick screening: there must be at least one carbon-carbon double bond somewhere.
    dbl_bond_smarts = Chem.MolFromSmarts("[C]=[C]")
    if not mol.HasSubstructMatch(dbl_bond_smarts):
        return False, "No carbon–carbon double bond (C=C) found in molecule"
    
    # Define SMARTS patterns for free acid and acyl ester.
    free_acid_smarts = Chem.MolFromSmarts("C(=O)[O;H1,O-]")  # free carboxylic acid (COOH or COO-)
    acyl_ester_smarts = Chem.MolFromSmarts("C(=O)O[C]")       # acyl ester motif: fatty acid in a lipid
    
    # Minimum length of the acyl chain (number of carbon atoms in the chain beyond the carbonyl)
    MIN_CHAIN_LENGTH = 6

    # This flag will be set if at least one qualifying acyl chain is found.
    chain_found = False
    reasons = []
    
    # Check free acid matches.
    free_matches = mol.GetSubstructMatches(free_acid_smarts)
    if free_matches:
        for match in free_matches:
            # In the free acid SMARTS "C(=O)[O;H1,O-]", match[0] is the carbonyl carbon.
            carboxyl_c = mol.GetAtomWithIdx(match[0])
            # For a free fatty acid, the acid group is terminal; the alpha carbon is any carbon neighbor
            # of the carboxyl carbon that is not the oxygen in the acid.
            alpha_c_found = False
            for nbr in carboxyl_c.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    alpha_c_found = True
                    start_atom = nbr
                    # Get the longest contiguous carbon chain from this starting point.
                    chain = longest_carbon_chain(start_atom, mol)
                    if len(chain) < MIN_CHAIN_LENGTH:
                        reasons.append(f"Free acid chain too short (length {len(chain)} < {MIN_CHAIN_LENGTH})")
                        continue
                    # Check that at least one bond in the chain is a double bond.
                    if not chain_has_double_bond(chain, mol):
                        reasons.append("Free acid chain does not contain a C=C double bond")
                        continue
                    # If this chain qualifies, we classify the molecule as an olefinic fatty acid.
                    return True, ("Contains a fatty acyl chain (via free acid) with sufficient length and a C=C double bond.")
            if not alpha_c_found:
                reasons.append("Free acid group found but no alpha carbon attached")
                
    # Check acyl ester matches (for cases such as phospholipids).
    ester_matches = mol.GetSubstructMatches(acyl_ester_smarts)
    if ester_matches:
        for match in ester_matches:
            # In the acyl ester SMARTS "C(=O)O[C]", match[0] is the carbonyl carbon,
            # match[1] is the ester oxygen, and match[2] is the acyl chain starting carbon.
            if len(match) < 3:
                continue
            start_atom = mol.GetAtomWithIdx(match[2])
            chain = longest_carbon_chain(start_atom, mol)
            if len(chain) < MIN_CHAIN_LENGTH:
                reasons.append(f"Acyl ester chain too short (length {len(chain)} < {MIN_CHAIN_LENGTH})")
                continue
            if not chain_has_double_bond(chain, mol):
                reasons.append("Acyl ester chain does not contain a C=C double bond")
                continue
            return True, ("Contains a fatty acyl chain (via acyl ester) with sufficient length and a C=C double bond.")
    
    # If we get here, no qualifying fatty acyl chain was found.
    if reasons:
        return False, " ".join(reasons)
    else:
        return False, "No free acid or acyl ester substructure (with long enough chain and unsaturation) found"

# Example usage (you can remove or comment these lines when integrating the function into a larger codebase).
if __name__ == "__main__":
    test_smiles_list = [
        "O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O",  # Avenoleic acid (should be True)
        "CC\\C=C/C\\C=C/CCC\\C=C\\C=C\\C=C/CCCC(O)=O",  # Icosapentaenoic acid (should be True)
        "CC\\C=C\\C(O)=O",  # trans-pent-2-enoic acid (should be False due to short chain)
        "OC[C@H](COP(OCC[N+](C)(C)C)(=O)[O-])OC(=O)CCCCCCC/C=C\\CCCCCC"  # a phospholipid fatty acyl chain (should be True)
    ]
    for s in test_smiles_list:
        result, reason = is_olefinic_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*50}")