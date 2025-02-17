"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: Phosphatidylglycerol 
Definition: A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.
Criteria applied in this implementation:
  1. Contains exactly one phosphorus atom.
  2. Among the three oxygen substituents on phosphorus (ignoring the double‐bonded oxygen),
     exactly one O leads to a branch that is “non‐acylated” (i.e. its immediately attached carbon 
     does not bear a carbonyl) – we assign this as the phosphoglycerol headgroup.
  3. The remaining branch (attached via a single P–O bond) should be a diacylglycerol.
     To check this we perform a graph search (avoiding the phosphorus) and count the number 
     of acyl ester groups (look for the substructure OC(=O)); there must be exactly two.
If any of these conditions fail, a reason is given.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    
    Our approach:
      1. Parse the molecule and check that it contains exactly one phosphorus atom.
      2. For that phosphorus atom, look at each oxygen connected by a single bond.
         (Skip the oxygen that is double-bonded to P.)
         Among these, we designate an oxygen as belonging to the headgroup branch if 
         the next atom (usually a carbon) is not “acylated” (i.e. is not directly bound 
         to a carbonyl group). All other branches (with a nearby C=O) are considered part 
         of the acyl (diacylglycerol) portion.
      3. For the acyl branch, we do a breadth‐first search (avoiding phosphorus) to 
         collect the branch atoms. Then we count all acyl esters (“OC(=O)”) that lie 
         entirely within that branch. A true PG should have exactly 2 such acyl ester groups.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is a phosphatidylglycerol, False otherwise.
      str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for exactly one phosphorus atom (atomic number 15)
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom, found {len(phosphorus_atoms)}"
    p_atom = phosphorus_atoms[0]
    
    # Get bonds of phosphorus. We only consider single bonds;
    # skip the oxygen that is double-bonded (P=O).
    single_O_neighbors = []
    for bond in p_atom.GetBonds():
        # Only consider bonds that are single
        if bond.GetBondType() == Chem.BondType.SINGLE:
            nbr = bond.GetOtherAtom(p_atom)
            if nbr.GetAtomicNum() == 8:  # oxygen
                single_O_neighbors.append(nbr)
    
    # We expect phosphorus to have exactly three substituents,
    # of which one should lead to the phosphoglycerol headgroup and
    # one should lead to the diacylglycerol (acyl) portion.
    if len(single_O_neighbors) < 2:
        return False, "Not enough oxygen substituents on phosphorus"
    
    headgroup_candidates = []
    acyl_candidates = []
    
    # For each oxygen neighbor (from the P-O single bonds), decide if
    # it is part of the headgroup or an acyl chain.
    # The idea: if the oxygen (let’s call it O_sub) is bonded to a carbon (C_sub)
    # that is directly connected to a carbonyl (C=O), then it is likely part of
    # an acyl chain. Otherwise, we treat it as the headgroup candidate.
    for o_atom in single_O_neighbors:
        # Get the neighbor (other than phosphorus)
        nbrs = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetIdx() != p_atom.GetIdx()]
        if not nbrs:
            continue
        # In typical lipids there will be only one other atom attached to this oxygen.
        c_atom = nbrs[0]
        # If c_atom is not carbon, skip further.
        if c_atom.GetAtomicNum() != 6:
            continue
        
        # Check if this carbon is directly involved in a carbonyl,
        # i.e. if any bond (except the one from o_atom) is a double bond to an oxygen.
        has_carbonyl = False
        for bond in c_atom.GetBonds():
            # skip the bond coming from o_atom
            if bond.GetBeginAtom().GetIdx() == o_atom.GetIdx() or bond.GetEndAtom().GetIdx() == o_atom.GetIdx():
                continue
            # look for a double bond to oxygen
            if (bond.GetBondType() == Chem.BondType.DOUBLE and
                (bond.GetBeginAtom().GetAtomicNum() == 8 or bond.GetEndAtom().GetAtomicNum() == 8)):
                has_carbonyl = True
                break
        
        if has_carbonyl:
            acyl_candidates.append(o_atom)
        else:
            headgroup_candidates.append(o_atom)
    
    if len(headgroup_candidates) != 1:
        return False, f"Phosphoglycerol headgroup branch not found exactly once (found {len(headgroup_candidates)} candidate(s))"
    if len(acyl_candidates) != 1:
        return False, f"Diacylglycerol (acyl) branch not found exactly once (found {len(acyl_candidates)} candidate(s))"
    
    # Now, for the acyl branch we would like to count the acyl ester groups (the –OC(=O)– fragments)
    # that belong to that branch. In a diacylglycerol we expect exactly 2.
    # To do this, we extract (via a breadth-first search) the set of atom indices that
    # belong to the acyl branch. We avoid traversing back to the phosphorus atom.
    def get_branch_atoms(start_atom, forbidden_idxs):
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
                # We allow traversal along any single bond.
                stack.append(nb)
        return branch

    # For the acyl branch, our "forbidden" set includes the phosphorus atom.
    forbidden = {p_atom.GetIdx()}
    acyl_branch_idxs = get_branch_atoms(acyl_candidates[0], forbidden)
    
    # Now count acyl ester groups in the entire molecule that reside in the acyl branch.
    # We define an ester group by the SMARTS: oxygen attached to a carbonyl: "OC(=O)"
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Only count those matches where the oxygen of the ester is part of the acyl branch.
    acyl_ester_count = 0
    for match in ester_matches:
        # match is a tuple of atom indices corresponding to (O, C, O)
        # we check if the first atom (the O) is in the acyl branch
        if match[0] in acyl_branch_idxs:
            acyl_ester_count += 1

    if acyl_ester_count != 2:
        return False, f"Expected 2 acyl ester groups in the diacylglycerol branch, found {acyl_ester_count}"
    
    # If all tests pass, we classify the molecule as a phosphatidylglycerol.
    return True, "Molecule contains one phosphorus atom with a unique phosphoglycerol headgroup branch and a diacylglycerol branch with 2 acyl esters."

# Example usage (testing with PG(8:0/8:0)):
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OC[C@@H](O)CO)(O)=O"
    result, reason = is_phosphatidylglycerol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)