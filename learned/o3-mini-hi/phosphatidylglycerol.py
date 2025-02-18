"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: Phosphatidylglycerol 
Definition: A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.
Criteria:
  1. The molecule must contain exactly one phosphorus atom.
  2. Among the oxygen neighbors attached with a single bond (ignoring dead-end oxygens
     and P=O bonds) exactly two branches are formed.
  3. One branch (the diacyl branch) must contain exactly two acyl ester fragments (defined
     by the pattern "OC(=O)" fully contained in the branch).
  4. The other branch (the phosphoglycerol headgroup) must contain no acyl ester fragments
     and should match a glycerol substructure pattern.
If any of these conditions are not met, a meaningful reason is returned.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.

    The algorithm:
      1. Parse the molecule and ensure exactly one phosphorus atom is present.
      2. For the phosphorus atom, retrieve oxygen neighbors that are connected via a single bond
         (ignoring double-bonded oxygens and dead ends).
      3. For each of the two candidate branches:
           a. Extract the branch as a sub-molecule (using a breadth-first search with the phosphorus
              as forbidden) so that SMARTS searches can be directly applied.
           b. Count the number acyl ester groups defined by the pattern "OC(=O)" that fall within the branch.
           c. If the branch contains 0 esters, check that it contains the glycerol headgroup motif
              via a SMARTS match (we use "OCC(O)CO" as a loose indicator).
           d. If the branch contains exactly 2 esters, mark it as the diacyl branch.
      4. If exactly one branch qualifies as the glycerol headgroup and exactly one as the diacyl
         branch, return True.
         Otherwise, return False along with an explanation.

    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a phosphatidylglycerol, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Ensure exactly one phosphorus atom (atomic number 15)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom, found {len(p_atoms)}"
    p_atom = p_atoms[0]
    
    # 2. Find oxygen neighbors attached by single bonds.
    candidate_oxygen_branches = []
    for bond in p_atom.GetBonds():
        # Only consider single bonds (skip P=O double bonds)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = bond.GetOtherAtom(p_atom)
        # Only consider oxygen atoms (atomic number 8)
        if neighbor.GetAtomicNum() != 8:
            continue
        # Ignore oxygen if it is a dead end (only connected to phosphorus)
        if neighbor.GetDegree() <= 1:
            continue
        candidate_oxygen_branches.append(neighbor)
        
    if len(candidate_oxygen_branches) != 2:
        return False, f"Expected 2 nontrivial oxygen branches from phosphorus, found {len(candidate_oxygen_branches)}"
    
    # Define helper function: breadth-first search to get the branch from a starting atom.
    def bfs_branch(start_atom, forbidden_idxs):
        branch = set()
        stack = [start_atom]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in branch:
                continue
            branch.add(atom.GetIdx())
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in forbidden_idxs:
                    continue
                stack.append(nb)
        return branch

    # We will not traverse back into the phosphorus atom.
    forbidden = {p_atom.GetIdx()}

    # Define ester SMARTS pattern: fragment "OC(=O)".
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    
    # Define glycerol headgroup SMARTS.
    # This pattern loosely captures a two-primary- hydroxyl glycerol connectivity:
    glycerol_smarts = "OCC(O)CO"  # note: this is a loose pattern and may be adjusted.
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    
    # Helper: Count acyl ester groups within a set of atom indices (branch)
    def count_acyl_esters(branch_idxs):
        count = 0
        # Get all matches of the ester pattern on the entire molecule.
        matches = mol.GetSubstructMatches(ester_pattern)
        for match in matches:
            # Count the ester only if the first atom (oxygen) is in the branch.
            if match[0] in branch_idxs:
                count += 1
        return count

    # Helper: Check if the branch submol has a glycerol headgroup fragment.
    def branch_is_glycerol(branch_idxs):
        try:
            branch_mol = Chem.PathToSubmol(mol, list(branch_idxs))
        except Exception as e:
            return False
        return branch_mol.HasSubstructMatch(glycerol_pattern)
    
    headgroup_branch_found = False
    acyl_branch_found = False
    headgroup_reason = ""
    acyl_reason = ""
    
    # Process each candidate branch.
    for o_atom in candidate_oxygen_branches:
        branch_idxs = bfs_branch(o_atom, forbidden)
        ester_count = count_acyl_esters(branch_idxs)
        if ester_count == 2:
            if acyl_branch_found:
                return False, "More than one branch with 2 acyl ester groups found"
            acyl_branch_found = True
            acyl_reason = f"Found diacyl branch with {ester_count} acyl ester group(s)."
        elif ester_count == 0:
            if branch_is_glycerol(branch_idxs):
                if headgroup_branch_found:
                    return False, "More than one branch appears to be a glycerol headgroup"
                headgroup_branch_found = True
                headgroup_reason = "Found glycerol headgroup branch (0 acyl ester groups with glycerol motif)."
            else:
                headgroup_reason += "A branch with 0 acyl esters did not show glycerol headgroup connectivity. "
        else:
            return False, f"Unexpected acyl ester count ({ester_count}) in one branch."
    
    if not headgroup_branch_found:
        return False, f"Phosphoglycerol headgroup branch not conclusively identified. {headgroup_reason}"
    if not acyl_branch_found:
        return False, f"Diacylglycerol branch (with exactly 2 acyl ester groups) not found. {acyl_reason}"
    
    return True, ("Molecule contains one phosphorus atom with a unique phosphoglycerol headgroup branch "
                  "and a diacylglycerol branch featuring exactly 2 acyl ester groups.")

# Example usage (you can test with one of the accepted PG SMILES):
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OC[C@@H](O)CO)(O)=O"  # PG(8:0/8:0)
    result, reason = is_phosphatidylglycerol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)