"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A cyclic fatty acid here is defined (for our purposes) as a molecule that
  (1) contains at least one ring,
  (2) does not contain amide bonds (to avoid peptides),
  (3) contains a terminal fatty acyl group – that is, a carbonyl functional group (acid or ester)
      with an oxygen (–OH or –OR) that is not attached to any nitrogen and is attached to exactly one
      acyclic, non‐aromatic carbon, and
  (4) starting from that attachment point a contiguous acyclic (and non‐aromatic) carbon “tail” (a linear chain)
      exists that has at least 5 carbons.
  
Note:
  By “terminal” we mean that the fatty acyl chain is attached by a single bond (and is unbranched)
  so that the tail can be “walked” linearly.
  
Examples of structures that should return True are listed in the long outcomes. If no valid fatty acyl group can be found
or if the tail is not long enough, then False is returned.
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is classified as a cyclic fatty acid based on its SMILES string.
    
    The molecule must have:
      - at least one ring present,
      - must not contain amide bonds (pattern: [NX3][CX3](=O)),
      - a fatty acyl group defined here as a carbonyl group [CX3](=O)[OX2H1,OX2] (i.e. acid or ester)
        that is not attached to any nitrogen, and from which there is exactly one non‐oxygen neighbor.
      - The non‐oxygen (typically a carbon) that is attached to the carbonyl must belong to a contiguous
        (unbranched) acyclic, non‐aromatic carbon “tail” of at least 5 atoms.
        
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria for a cyclic fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check if at least one ring exists.
    if not mol.GetRingInfo().AtomRings():
        return False, "No ring found in the structure"
    
    # 2. Reject molecules with amide bonds (avoid peptides).
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds (likely a peptide)"
    
    # 3. Identify fatty acyl group candidates.
    # SMARTS: [CX3](=O)[OX2H1,OX2] matches a carbonyl (acid or ester)
    fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX2]")
    valid_chain_lengths = []  # will collect chain lengths for valid fatty acyl groups
    for match in mol.GetSubstructMatches(fatty_acyl_pattern):
        # In the match, match[0] is the carbonyl carbon; match[1] is the oxygen attached.
        carbonyl = mol.GetAtomWithIdx(match[0])
        # Skip if the carbonyl is attached to any nitrogen (would be an amide)
        if any(nbr.GetAtomicNum() == 7 for nbr in carbonyl.GetNeighbors()):
            continue
        
        # Find non-oxygen neighbors of the carbonyl (should be exactly one for a terminal fatty acyl group).
        non_oxygen_neighbors = [nbr for nbr in carbonyl.GetNeighbors() if nbr.GetAtomicNum() != 8]
        if len(non_oxygen_neighbors) != 1:
            continue  # not a terminal acyl group
        start_atom = non_oxygen_neighbors[0]
        # We require that the start_atom is a carbon and is acyclic and non‐aromatic.
        if start_atom.GetAtomicNum() != 6 or start_atom.IsInRing() or start_atom.GetIsAromatic():
            continue

        # 4. Build the set of candidate acyclic carbons (non-ring, non‐aromatic carbons).
        acyclic_atoms = {atom.GetIdx() for atom in mol.GetAtoms()
                         if atom.GetAtomicNum() == 6 and (not atom.IsInRing()) and (not atom.GetIsAromatic())}
        if start_atom.GetIdx() not in acyclic_atoms:
            continue

        # For a terminal fatty acyl tail, we want the chain to be linear from the attachment point.
        # In an acyclic (non‐aromatic) carbon subgraph, build an adjacency list.
        neighbors_dict = {idx: [] for idx in acyclic_atoms}
        for idx in acyclic_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in acyclic_atoms:
                    neighbors_dict[idx].append(nbr.GetIdx())
        
        # In a terminal chain, the first carbon (start_atom) should have degree 1 in the acyclic subgraph.
        if len(neighbors_dict[start_atom.GetIdx()]) != 1:
            # If it is branched, we skip it (for our strict definition).
            continue
        
        # Now walk linearly from the start_atom.
        chain_length = 1  # count start_atom as first in chain
        current = start_atom.GetIdx()
        previous = None
        while True:
            # Get neighbors (in the acyclic graph) excluding the one we came from.
            nbrs = [n for n in neighbors_dict[current] if n != previous]
            if len(nbrs) == 1:
                chain_length += 1
                previous = current
                current = nbrs[0]
            else:
                break
        
        valid_chain_lengths.append(chain_length)
    
    if not valid_chain_lengths:
        return False, "No valid terminal fatty acyl group found (acid/ester with a single attachment to an acyclic carbon)"
    
    # For our purpose, we require that at least one acyl tail is of length >= 5.
    max_chain = max(valid_chain_lengths)
    if max_chain < 5:
        return False, f"No sufficiently long contiguous acyclic carbon chain (>=5 carbons) found; found chain length = {max_chain}"
    
    return True, f"Contains a terminal fatty acyl group, a contiguous acyclic carbon chain of length {max_chain}, and at least one ring in the structure"

# For testing:
if __name__ == "__main__":
    # Some examples taken from the outcomes:
    examples = [
        # True positives:
        ("O1C(CCCCCCCCCCC(OC[C@@H](OC(=O)CCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C", "DG(11D5/11M5/0:0)"),
        ("OC(=O)CCCCCCNC1c2ccccc2CCc2ccccc12", "amineptine"),
        ("C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O", "(5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid"),
        ("C(CCC)C/C=C\\C[C@H]1[C@@H](C/C=C\\C/C=C\\CCCC(O)=O)O1", "(11S,12R)-EET"),
        ("[H]C(\\C=C/CCCCCCCC(O)=O)=C1O[C@H]1C\\C=C/CC", "(9Z,13S,15Z)-12,13-epoxyoctadeca-9,11,15-trienoic acid"),
        # False negatives examples (which our altered algorithm may still miss):
        ("O=C1OCC(=C1C)CCCC(=O)O", "Ficifuranone A"),
        ("O1C(CCCCC(O)=O)=C(C=C1CCCCC)", "3-Methyl-5-pentyl-2-furanpentanoic acid"),
        ("O1C(=CC=C1C(=O)C#C/C=C/C)/C=C\\C(O)=O", "Wyeronic acid"),
        ("OC(=O)CCCC[C@@H]1CCSS1", "(R)-lipoic acid")
    ]
    for smi, name in examples:
        result, reason = is_cyclic_fatty_acid(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")