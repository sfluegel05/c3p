"""
Classifies: CHEBI:87657 octanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Octanoate Ester
Definition: Any fatty acid ester in which the carboxylic acid component 
is octanoic acid (caprylic acid). In other words, every ester group in the 
molecule must have an acyl chain that is exactly CH3–(CH2)6–C(=O)O.
This implementation first locates ester groups (using a general SMARTS for C(=O)O)
and then, for each ester, “walks” from the carbonyl carbon along the acyl chain 
to check that it is a linear chain with exactly 8 carbons (the carbonyl plus 7 alkyl carbons) 
with no branching.
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if every ester group in the molecule is derived from octanoic acid.
    
    For each ester group (identified by a C(=O)O pattern), the algorithm finds the
    acyl chain connected to the carbonyl carbon and traverses it. The ester qualifies as 
    an octanoate ester if the acyl chain is linear and exactly eight carbons long 
    (i.e. the carbonyl plus 7 additional carbon atoms, corresponding to CH3–(CH2)6–C(=O)).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if every ester group is derived from octanoic acid, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to find a general ester group: a carbonyl carbon connected to an oxygen.
    ester_pattern = Chem.MolFromSmarts("[C:1](=O)[O:2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern, uniquify=True)
    if not ester_matches:
        return False, "No ester group found"
    
    for match in ester_matches:
        carbonyl_idx, oxy_idx = match
        # For each ester, check that the acyl chain (from the carbonyl carbon) is exactly octanoate.
        if not _check_octanoate_chain(mol, carbonyl_idx, oxy_idx):
            return False, "Found an ester group whose acyl chain is not derived from octanoic acid"
    
    return True, "All ester groups are derived from octanoic acid (octanoate ester)."


def _check_octanoate_chain(mol, carbonyl_idx: int, ester_oxy_idx: int) -> bool:
    """
    Given an ester group identified by its carbonyl carbon (carbonyl_idx)
    and the oxygen forming the ester bond (ester_oxy_idx), verify that the acyl chain 
    attached to the carbonyl (i.e. on the acid side) is exactly octanoate:
        CH3–CH2–CH2–CH2–CH2–CH2–CH2–C(=O)O

    The algorithm:
      1. From the carbonyl carbon, get the neighbor that is a carbon
         (ignoring the ester oxygen).
      2. "Walk" along the acyl chain while enforcing linearity (no branching).
      3. Count the total number of carbons in the acyl chain INCLUDING the carbonyl.
         For octanoate this must equal 8.
      4. Also ensure that the terminal atom is a CH3 (only one carbon–carbon neighbor).

    Returns:
        bool: True if the acyl chain is exactly octanoic acid derived; False otherwise.
    """
    carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the acyl neighbor on the acid side.
    acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() 
                      if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != ester_oxy_idx]
    # It should have exactly one acyl neighbor.
    if len(acyl_neighbors) != 1:
        return False
    acyl_start = acyl_neighbors[0]
    
    # Traverse the acyl chain starting from acyl_start.
    chain_atoms = _traverse_linear_chain(mol, start_atom=acyl_start, parent=carbonyl)
    if chain_atoms is None:
        # Branching detected or other issue.
        return False
    # Total number of carbons from the carbonyl plus the chain atoms.
    total_carbons = 1 + len(chain_atoms)
    # Octanoate must have exactly 8 carbons.
    if total_carbons != 8:
        return False
    # Check that the terminal atom (last in chain_atoms) is a methyl group: 
    # it should have only one carbon neighbor.
    terminal = chain_atoms[-1]
    carbon_neighbors = [nbr for nbr in terminal.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False
    return True


def _traverse_linear_chain(mol, start_atom, parent):
    """
    Traverses a linear carbon chain starting from start_atom.
    'parent' is the atom from which we came (to avoid backtracking).
    
    Returns:
        list: a list of atoms encountered (in order) along the chain.
              If branching is detected (more than one carbon neighbor excluding the parent),
              returns None.
    """
    chain = []
    current = start_atom
    prev = parent
    chain.append(current)
    while True:
        # Look for carbon neighbors of 'current' excluding the one we came from.
        nbrs = [nbr for nbr in current.GetNeighbors() 
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev.GetIdx()]
        if len(nbrs) == 0:
            # Reached terminal carbon.
            break
        if len(nbrs) > 1:
            # Branching detected: not a simple linear chain.
            return None
        # Advance along the chain.
        prev, current = current, nbrs[0]
        chain.append(current)
    return chain


# The code below is for testing purposes only.
if __name__ == '__main__':
    test_smiles = [
        # True positives (should return True):
        "CCCCCCCC(=O)OC[C@H](O)CO",  # 3-octanoyl-sn-glycerol
        "CCCCCCCC(=O)OC",            # methyl octanoate
        "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCC",  # 1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate)
        "CCCCCCCC(=O)OC[C@H](O)CO",  # 1-octanoyl-sn-glycerol (similar to above)
        "CCCCCCCC(=O)OCC",           # ethyl octanoate
        "CCCCCCCC(=O)OCCC",          # propyl octanoate
        # False positives / negatives (should return False):
        "O(C(=O)CCCCCCC)C(C)C",      # isopropyl octanoate (WRONGLY CLASSIFIED in previous attempt)
        "O(CCCCCCCCCCCC)C(=O)CCCCCCC",# dodecyl octanoate (acyl chain too long)
        "CCCCCCCC(=O)O[C@H](COP(O)(=O)OP(O)(O)=O)OC(=O)CCCCCCCC",  # 1,2-dioctanoyl-sn-glycerol 3-diphosphate is expected True,
                                                      # but other false negatives should occur if an ester group isn’t pure octanoate.
    ]
    for sm in test_smiles:
        flag, reason = is_octanoate_ester(sm)
        print(f"SMILES: {sm}\nClassification: {flag}\nReason: {reason}\n")