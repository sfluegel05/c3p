"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: Phosphatidylglycerol 
Definition: A glycerophosphoglycerol that is glycerol in which one hydrogen from a primary hydroxy group has been replaced by a phosphatidyl group.
Criteria:
  1. The molecule must contain exactly one phosphorus atom.
  2. Of the oxygen atoms connected to phosphorus by single bonds (ignoring P=O bonds and dead‐end oxygens), exactly two branches should be present.
  3. One branch (the diacyl branch) must contain exactly 2 acyl ester groups (the “OC(=O)” fragment).
  4. The other branch (the phosphoglycerol headgroup) must contain 0 acyl ester groups and match a glycerol-like substructure.
If any condition is not met, a meaningful reason is returned.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines whether a molecule is a phosphatidylglycerol based on its SMILES string.

    The algorithm:
      1. Parse the molecule and ensure exactly one phosphorus atom is present.
      2. For the phosphorus atom, examine its oxygen neighbors connected via single bonds
         (ignoring double-bonded oxygens and "dead-end" oxygens having no further neighbors apart from P).
      3. For each oxygen branch, perform a breadth-first search (with a depth limit)
         to extract a local subgraph of the branch.
      4. Count acyl ester groups (using the SMARTS "OC(=O)") within that branch.
         If a branch contains exactly 2 esters, it is marked as the diacyl (acyl) branch.
      5. If a branch contains 0 esters and its sub-molecule is acyclic and matches a loose glycerol 
         headgroup pattern (here “OCC(O)CO”), then it is taken as the phosphoglycerol headgroup branch.
      6. If exactly one branch qualifies for each role, the molecule is classified as a phosphatidylglycerol.

    Args:
        smiles (str): A SMILES string.

    Returns:
        bool: True if the molecule is classified as phosphatidylglycerol; otherwise False.
        str: A message indicating the reason for the decision.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check that there is exactly one phosphorus atom (atomic number 15)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom, found {len(p_atoms)}"
    p_atom = p_atoms[0]
    
    # 2. Find candidate oxygen branches: oxygen atoms connected to phosphorus via a single bond,
    # ignoring oxygens that are dead ends (degree<=1) and ignoring bonds that are not single.
    candidate_oxygens = []
    for bond in p_atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = bond.GetOtherAtom(p_atom)
        if neighbor.GetAtomicNum() != 8:
            continue
        if neighbor.GetDegree() <= 1:
            continue
        candidate_oxygens.append(neighbor)
    
    if len(candidate_oxygens) != 2:
        return False, f"Expected 2 nontrivial oxygen branches from phosphorus, found {len(candidate_oxygens)}"
    
    # Helper: breadth-first search on the molecular graph with a depth limit.
    def bfs_branch(start_atom, forbidden_idxs, max_depth=6):
        visited = set()
        queue = [(start_atom, 0)]
        while queue:
            atom, depth = queue.pop(0)
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if depth < max_depth:
                for nb in atom.GetNeighbors():
                    if nb.GetIdx() in forbidden_idxs:
                        continue
                    if nb.GetIdx() not in visited:
                        queue.append((nb, depth + 1))
        return visited
    
    # We will not traverse back into the phosphorus atom.
    forbidden = {p_atom.GetIdx()}
    
    # Define SMARTS patterns.
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    
    # A loose glycerol headgroup indicator.
    glycerol_smarts = "OCC(O)CO"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    
    # Helper: Count the acyl ester groups within the branch (only count if the oxygen in "OC(=O)" is in the branch).
    def count_acyl_esters(branch_idxs):
        count = 0
        matches = mol.GetSubstructMatches(ester_pattern)
        for match in matches:
            # We expect the match to start with the oxygen atom of the ester fragment.
            if match[0] in branch_idxs:
                count += 1
        return count
    
    # Helper: Check if the branch sub-molecule matches a glycerol headgroup.
    # We also require that the branch be acyclic to avoid inadvertently capturing larger ring systems.
    def branch_matches_glycerol(branch_idxs):
        try:
            branch_mol = Chem.PathToSubmol(mol, list(branch_idxs))
        except Exception:
            return False
        if branch_mol.GetRingInfo().NumRings() > 0:
            return False
        return branch_mol.HasSubstructMatch(glycerol_pattern)
    
    headgroup_found = False
    acyl_found = False
    headgroup_reason = ""
    acyl_reason = ""
    
    # Process each candidate oxygen branch.
    for o_atom in candidate_oxygens:
        branch_idxs = bfs_branch(o_atom, forbidden, max_depth=6)
        ester_count = count_acyl_esters(branch_idxs)
        
        if ester_count == 2:
            # Diacyl (acyl) branch identified.
            if acyl_found:
                return False, "More than one branch with 2 acyl ester groups found"
            acyl_found = True
            acyl_reason = f"Found diacyl branch with {ester_count} acyl ester group(s)."
        elif ester_count == 0:
            # Check if it fulfills the glycerol headgroup criteria.
            if branch_matches_glycerol(branch_idxs):
                if headgroup_found:
                    return False, "More than one branch appears to be a glycerol headgroup"
                headgroup_found = True
                headgroup_reason = "Found glycerol headgroup branch (0 acyl esters; glycerol motif confirmed)."
            else:
                headgroup_reason += "A branch with 0 acyl esters did not match the glycerol motif. "
        else:
            return False, f"Unexpected acyl ester count ({ester_count}) in one branch."
    
    if not headgroup_found:
         return False, "Phosphoglycerol headgroup branch not conclusively identified. " + headgroup_reason
    if not acyl_found:
         return False, "Diacylglycerol branch (with exactly 2 acyl ester groups) not found. " + acyl_reason
    
    return True, ("Molecule contains one phosphorus atom with a unique phosphoglycerol headgroup branch "
                   "and a diacylglycerol branch featuring exactly 2 acyl ester groups.")

# Example usage
if __name__ == "__main__":
    # Test with one of the valid PG examples.
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OC[C@@H](O)CO)(O)=O"  # PG(8:0/8:0)
    result, reason = is_phosphatidylglycerol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)