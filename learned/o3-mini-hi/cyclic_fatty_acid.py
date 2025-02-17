"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A cyclic fatty acid here is expected to have a fatty acyl functional group (either a terminal carboxylic acid
or ester function not attached on both sides to carbon), a contiguous acyclic (non‐ring, non‐aromatic) carbon chain of length >= 5,
and at least one ring somewhere in the structure.
Additionally, molecules containing any amide bond (common in peptides) are rejected.
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is classified as a cyclic fatty acid based on its SMILES string.
    
    The molecule must have:
      - at least one ring (of any size),
      - a terminal fatty acyl group. In our filter, a valid fatty acyl group is a carbonyl
        [CX3](=O)[OX2H1,OX2] that is not attached to any nitrogen and is attached to exactly one non‐oxygen carbon.
      - from the carbon attached to the fatty acyl group, a contiguous chain of 5 or more acyclic,
        non‐aromatic carbon atoms (this chain does not include the carbonyl carbon itself).
      - No amide bonds (identified via the pattern [NX3][CX3](=O)) – this helps avoid peptides.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a cyclic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for presence of at least one ring.
    if not mol.GetRingInfo().AtomRings():
        return False, "No ring found in the structure"
    
    # 2. Check for amide bonds (common in peptides) and if found, reject.
    #    We search for a pattern where a nitrogen is directly bonded to a carbonyl carbon.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if amide_pattern is None:
        return False, "Failed to generate amide SMARTS pattern"
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds (likely a peptide)"
    
    # 3. Identify a valid fatty acyl functional group.
    # Accepts both acids and esters by allowing the second oxygen to be hydroxyl (acid) or bound to carbon (ester).
    # We use the SMARTS pattern for a carbonyl with an oxygen: [CX3](=O)[OX2H1,OX2]
    fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX2]")
    if fatty_acyl_pattern is None:
        return False, "Failed to generate fatty acyl SMARTS pattern"
    
    valid_start_atoms = []
    for match in mol.GetSubstructMatches(fatty_acyl_pattern):
        # match[0] is the carbonyl carbon; match[1] is the oxygen bound to it.
        carbonyl = mol.GetAtomWithIdx(match[0])
        # Check that carbonyl is not attached to any nitrogen (i.e. not an amide)
        if any(nbr.GetAtomicNum() == 7 for nbr in carbonyl.GetNeighbors()):
            continue
        # Now collect non-oxygen neighbors (i.e. carbons) that are candidates for the fatty acyl chain.
        chain_neighbors = []
        for nbr in carbonyl.GetNeighbors():
            # Exclude neighbors that are oxygen (which are already part of the carbonyl or –OH)
            if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()) and (not nbr.GetIsAromatic()):
                chain_neighbors.append(nbr.GetIdx())
        if len(chain_neighbors) == 1:
            valid_start_atoms.extend(chain_neighbors)
    
    if not valid_start_atoms:
        return False, "No terminal fatty acyl group found (valid acid or ester functionality not detected)"
    
    # 4. Compute the longest contiguous acyclic aliphatic carbon chain
    # We'll restrict DFS to carbons which are:
    #   - AtomicNum==6
    #   - Not in any ring
    #   - Not aromatic
    acyclic_carbons = {atom.GetIdx() for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 and (not atom.IsInRing()) and (not atom.GetIsAromatic())}
    
    # Build an adjacency list (graph) for the acyclic carbon subgraph.
    neighbors = {idx: [] for idx in acyclic_carbons}
    for idx in acyclic_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in acyclic_carbons:
                neighbors[idx].append(nbr.GetIdx())
    
    # DFS to compute the longest contiguous chain starting from a particular start atom.
    def dfs(current, visited):
        max_length = 1
        for nbr in neighbors[current]:
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    longest_chain = 0
    # We require that at least one chain starting at a valid fatty acyl attachment point is >= 5.
    for start in valid_start_atoms:
        if start not in acyclic_carbons:
            continue
        chain_length = dfs(start, {start})
        if chain_length > longest_chain:
            longest_chain = chain_length

    if longest_chain < 5:
        return False, f"No sufficiently long contiguous aliphatic chain (>=5 carbons) found; found chain length = {longest_chain}"
    
    return True, f"Contains a terminal fatty acyl group, a contiguous acyclic carbon chain of length {longest_chain}, and at least one ring in the structure"

# For testing:
if __name__ == "__main__":
    # Running a few examples from the provided outcomes.
    test_examples = [
        # Expected true positives:
        ("OC(=O)CCCCCCNC1c2ccccc2CCc2ccccc12", "amineptine"),
        ("C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O", "(5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid"),
        ("O1C(CCCCC(O)=O)=C(C=C1CCCCC)C", "3-Methyl-5-pentyl-2-furanpentanoic acid"),
        ("OC(=O)CCCC[C@@H]1CCSS1", "(R)-lipoic acid"),
        # One false negative (should be cyclic fatty acid but was missed earlier):
        ("O1C(CCCCCCCCCCC(OC[C@@H](OC(=O)CCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C)C", "DG(11D5/11M5/0:0)")
    ]
    
    for sm, name in test_examples:
        res, reason = is_cyclic_fatty_acid(sm)
        print(f"Name: {name}\nSMILES: {sm}\nResult: {res}\nReason: {reason}\n{'-'*60}")