"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated Fatty Acid
Definition: A fatty acid is defined here as a molecule having a terminal carboxylic acid group 
attached to a single, linear, non‐aromatic aliphatic chain that contains more than one 
non‐aromatic C=C double bond. Additional checks attempt to ensure that the chain is unbranched 
and that no extra functional groups (for example, amides) are present along the chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a polyunsaturated fatty acid (PUFA) based on:
      1. Presence of a terminal carboxylic acid group, where the acid carbon is attached to exactly one carbon.
      2. A linear (“walkable”), unbranched, acyclic, non‐aromatic aliphatic chain traced from that single neighbor.
         (Any branching in the carbon chain will cause rejection.)
      3. The absence of any carbon–carbon triple bonds along the traced chain.
      4. The presence of at least two non‐aromatic C=C double bonds (counted along the chain).
      5. A couple of extra heuristics: minimal chain length and no extra (non‐Heteroatom) substituents and no amide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as a PUFA, otherwise False.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 0. Reject molecules with amide bonds (outside the carboxylic acid). We define an amide
    # SMARTS that should not be present (excluding the terminal -COOH itself).
    amide_smarts = Chem.MolFromSmarts("[NX3;!$(NC(=O)[O;H])][CX3](=O)[#6]")  # any N-C(=O)-C (roughly)
    if mol.HasSubstructMatch(amide_smarts):
        return False, "Molecule contains an amide bond which is not allowed in a typical fatty acid"
    
    # 1. Look for carboxylic acid group
    # We use SMARTS matching a COOH: carbon with a double-bonded O and an OH.
    acid_smarts = "[CX3](=O)[O;H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No terminal carboxylic acid group found; not a fatty acid"
    
    # 2. Among acid groups, pick one where the acid carbon (first atom in match) is attached to exactly one carbon.
    terminal_acid_found = False
    acid_atom = None
    chain_start = None
    for match in acid_matches:
        # We assume match[0] is the acid carbon (sp2 C of the COOH)
        acid_c = mol.GetAtomWithIdx(match[0])
        # Find neighboring carbon atoms (atomic num 6)
        carbon_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acid_atom = acid_c
            chain_start = carbon_neighbors[0]
            break
    if not terminal_acid_found or chain_start is None:
        return False, "Carboxylic acid group is not terminal; not a typical fatty acid"
    
    # 3. Trace the acyl chain beginning from chain_start.
    # The rule is: at each step, only one non‐aromatic, non‐ring carbon (atomic num 6) that is not the atom we came from is accepted.
    chain_atom_idxs = []  # indices in order along the chain
    current_atom = chain_start
    previous_atom = acid_atom
    chain_atom_idxs.append(current_atom.GetIdx())
    
    while True:
        # Look for neighbor carbons that are (a) not the atom we came from, (b) non‐aromatic, and (c) not in any ring.
        candidates = []
        for nbr in current_atom.GetNeighbors():
            if nbr.GetIdx() == previous_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            candidates.append(nbr)
        if len(candidates) == 0:
            # End of chain reached.
            break
        if len(candidates) > 1:
            return False, "Acyl chain is branched; expected a single linear chain for a typical fatty acid"
        # Proceed along the only available neighbor.
        next_atom = candidates[0]
        chain_atom_idxs.append(next_atom.GetIdx())
        previous_atom = current_atom
        current_atom = next_atom

    chain_length = len(chain_atom_idxs)
    MIN_CHAIN_LENGTH = 8  # at least 8 carbons in the chain (this is a heuristic)
    if chain_length < MIN_CHAIN_LENGTH:
        return False, f"Fatty acid chain length only {chain_length} carbons; too short to be typical"
    
    # 4. Verify that the chain is truly unbranched.
    # For each carbon in the chain (except the first one which is attached to the acid),
    # check that no extra carbon-neighbor (other than the one coming from the chain) is present.
    # (Extra carbon substituents would indicate branching.)
    chain_atom_set = set(chain_atom_idxs)
    for idx in chain_atom_idxs:
        atom = mol.GetAtomWithIdx(idx)
        # Consider all carbon neighbors that are non‐aromatic and not in a ring.
        ext_cntrs = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                # Allow neighbor if it is either the preceding or following atom in the chain.
                if nbr.GetIdx() in chain_atom_set:
                    continue
                ext_cntrs.append(nbr)
        if ext_cntrs:
            return False, "Acyl chain is substituted (has extra carbon substituents); expected a single linear chain"
    
    # 5. Now count the non‐aromatic double bonds and check for any triple bonds along the (traced) chain.
    double_bond_count = 0
    for i in range(len(chain_atom_idxs)-1):
        a1 = chain_atom_idxs[i]
        a2 = chain_atom_idxs[i+1]
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is None:
            continue
        # Reject if the bond is a triple bond.
        if bond.GetBondType() == Chem.BondType.TRIPLE:
            return False, "Chain contains a carbon–carbon triple bond which is not allowed in typical fatty acids"
        # Count double bonds that are not aromatic.
        if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
            double_bond_count += 1

    if double_bond_count < 2:
        return False, f"Found {double_bond_count} non-aromatic C=C bond(s) in chain; need at least 2 for polyunsaturation"
    
    # All checks passed: provide a reason message that includes chain length and double bond count.
    reason = (f"Contains a terminal carboxylic acid group attached to a linear chain of {chain_length} carbons "
              f"with {double_bond_count} non-aromatic double bond(s); qualifies as a polyunsaturated fatty acid")
    return True, reason


# Example usage:
if __name__ == "__main__":
    # Test with one example; for instance, resolvin D6.
    test_smiles = "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O"
    result, explanation = is_polyunsaturated_fatty_acid(test_smiles)
    print(result, explanation)