"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: Any branched‐chain fatty acid containing methyl branches only.
The function is_methyl_branched_fatty_acid takes a SMILES string as input and returns a boolean plus a reason.
The approach:
  1. Parse molecule, require acyclic and presence of a carboxylic acid group.
  2. Build a carbon–only graph from which we enumerate all simple paths beginning at the acid carbon.
  3. For each candidate path (of at least 5 carbons) we check:
       – that any additional carbon neighbors (i.e. branches) are terminal – thus are CH3 groups.
       – that any non‐carbon neighbor attached to a backbone carbon is acceptable only when:
            • on the acid carbon: any oxygen is allowed,
            • on the terminal (non–acid) carbon: an –OH (oxygen with one hydrogen) is allowed,
            • on any interior backbone carbon: no heteroatom is allowed.
       – additionally, if a backbone carbon has a nitrogen substituent it is only allowed if the N has no hydrogens
            (as in a quaternary ammonium, e.g. in a betaine) – this helps to reject amino acid cases.
 4. Among acceptable candidates we choose the one with maximal main–chain length.
If none passes the tests then an error message is returned.
"""

from rdkit import Chem
from collections import deque

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl–branched fatty acid (branched–chain fatty acid with methyl branches only)
    based on its SMILES string.
    
    Args:
      smiles (str): SMILES string.
    
    Returns:
      bool: True if classified as methyl–branched fatty acid.
      str: Explanation/reason for classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Only allow acyclic molecules (fatty acids are expected to be open–chain).
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a straight–chain fatty acid"
    
    # Look for a carboxylic acid group.
    # The SMARTS below matches protonated or deprotonated –COOH groups.
    acid_smarts = "C(=O)[O;H,-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    # We choose the acid carbon from the first match.
    acid_carbon = acid_matches[0][0]
    
    # For checking hydrogen counts we add explicit H to the molecule.
    mol_h = Chem.AddHs(mol)
    
    # Build a carbon-only graph from the explicit–H molecule.
    carbon_indices = {atom.GetIdx() for atom in mol_h.GetAtoms() if atom.GetAtomicNum() == 6}
    if acid_carbon not in carbon_indices:
        return False, "Carboxyl carbon not found among carbons"
    carbon_adj = {idx: set() for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol_h.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_adj[idx].add(nbr.GetIdx())
    
    # Enumerate all simple paths in the carbon graph starting from the acid carbon.
    all_paths = []
    def dfs(current, path):
        extended = False
        for nbr in carbon_adj[current]:
            if nbr in path:
                continue
            dfs(nbr, path + [nbr])
            extended = True
        if not extended:
            all_paths.append(path)
    dfs(acid_carbon, [acid_carbon])
    
    if not all_paths:
        return False, "No carbon chain found"
    
    # Helper function to decide if a heteroatom attached to a backbone carbon is acceptable.
    def is_allowed_substituent(backbone_idx, nbr_atom, terminal=False):
        sym = nbr_atom.GetSymbol()
        # For the acid carbon, allow any oxygen (the acid group normally has two oxygens).
        if backbone_idx == acid_carbon:
            if sym == 'O':
                return True
            return False
        # For the terminal backbone carbon (other than acid carbon), allow an -OH group.
        if terminal:
            if sym == 'O':
                # Expect it to be hydroxyl (one H neighbor).
                num_H = sum(1 for n in nbr_atom.GetNeighbors() if n.GetSymbol()=='H')
                if num_H == 1:
                    return True
            # Also allow a quaternary nitrogen group here (if it has no H’s)
            if sym == 'N' and nbr_atom.GetTotalNumHs() == 0:
                return True
            return False
        # For interior backbone carbons, do not allow any non-carbon substituents
        # except (optionally) a nitrogen that is quaternary.
        if sym == 'N' and nbr_atom.GetTotalNumHs() == 0:
            return True
        return False
    
    # Check each candidate main chain (a simple path from the acid carbon).
    # We require that the chain has at least 5 carbons.
    # For each backbone carbon (each atom in the chain) we check:
    #   (a) all attached atoms that are NOT carbons are allowed (using our simple rules);
    #   (b) all carbon-atom attachments not in the chain are valid methyl branches.
    # A valid methyl branch (as seen in the carbon graph) must be terminal (degree 1) and must, after H-addition,
    # have exactly 3 hydrogens.
    valid_candidates = []
    def check_chain(chain):
        if len(chain) < 5:
            return False, "Main carbon chain is too short to be a fatty acid", None
        branch_count = 0
        chain_set = set(chain)
        # For each backbone carbon in the candidate chain…
        for i, idx in enumerate(chain):
            atom = mol_h.GetAtomWithIdx(idx)
            # Check non-carbon neighbors (heteroatoms).
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    continue
                # Only consider atoms NOT already part of a valid acid group.
                # For acid carbon, any oxygen is allowed.
                # For terminal backbone atoms (last in chain, and not the acid carbon) allow -OH or a quaternary N.
                if idx == chain[0]:
                    if not (nbr.GetSymbol() == 'O'):
                        return False, f"Acid carbon has disallowed substituent {nbr.GetSymbol()}", None
                elif idx == chain[-1]:
                    if not is_allowed_substituent(idx, nbr, terminal=True):
                        return False, f"Terminal backbone carbon has disallowed substituent {nbr.GetSymbol()}", None
                else:
                    if not is_allowed_substituent(idx, nbr, terminal=False):
                        return False, f"Interior backbone carbon has disallowed substituent {nbr.GetSymbol()}", None
            # Now check carbon branches: take neighbors in the carbon graph that are not in the backbone.
            for nbr_idx in carbon_adj[idx]:
                if nbr_idx in chain_set:
                    continue
                # The branch (nbr_idx) must be terminal in the carbon graph.
                if len(carbon_adj[nbr_idx]) != 1:
                    return False, f"Branch at backbone atom {idx} is not a terminal carbon", None
                branch_atom = mol_h.GetAtomWithIdx(nbr_idx)
                # Check that the branch carbon is “methyl‐like”: exactly 3 hydrogens.
                if branch_atom.GetTotalNumHs() != 3:
                    return False, f"Branch carbon at index {nbr_idx} is not CH3", None
                branch_count += 1
        return True, "", branch_count

    candidate_chains = []
    for chain in all_paths:
        valid, err, branches = check_chain(chain)
        if valid:
            candidate_chains.append((chain, branches))
    if not candidate_chains:
        # If none pass, report error from the longest chain candidate.
        longest = max(all_paths, key=len)
        if len(longest) < 5:
            return False, "Main carbon chain is too short to be a fatty acid"
        return False, "No candidate main chain found that yields only methyl branches"
    # Among candidates, choose the chain with maximal length.
    best_chain, best_branch_count = max(candidate_chains, key=lambda x: len(x[0]))
    msg = (f"CORRECT Methyl-branched fatty acid with {best_branch_count} methyl branch(es) "
           f"on a main chain of length {len(best_chain)}")
    return True, msg

# Example usage:
if __name__ == '__main__':
    test_smiles_list = [
        "CC(C)CCCCCCCCCCCCCCCCC(O)=O",           # 18-methylnonadecanoic acid - expect True
        "CC(CCCCCCC/C=C/C(=O)O)C",                # (E)-11-methyldodec-2-enoic acid - expect True
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",      # 24-methylpentacosanoic acid - expect True
        "OC(=O)CC(C)C=C",                        # 3-methyl-4-pentenoic acid - expect True
        "CC(C)CCCCCCCCC(O)=O",                    # 10-methylundecanoic acid - expect True
        "OC(=O)CCC(C)=C",                        # 4-methyl-4-pentenoic acid - expect True
        "C(CCCCCCC(CC)C)CCCCC(O)=O",              # 13-methylpentadecanoic acid - expect True
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methylnonacosanoic acid - expect True
        "CC(C)CCCCCCCCCCCCCCCCCCC(O)=O",          # 20-methylhenicosanoic acid - expect True
        "CC(CO)CCCCCCCCCCCCCC(O)=O",              # omega-hydroxy-15-methylpalmitic acid - expect True
        "OC(CCCCCCCCCCCCCCCC(C)C)=O",             # 17-methyloctadecanoic acid - expect True
        "CCC(C)CCCCCCCCCCCCCCCCCCCCC(O)=O",        # 22-methyltetracosanoic acid - expect True
        "[O-]C(=O)CCC([N+](C)(C)C)C",              # 4-aminovaleric acid betaine - expect True
        "OC(=O)/C=C(\\CC)/C",                     # 3-methyl-2Z-pentenoic acid - expect True
        "CCC(C)CC(O)=O",                         # 3-methylvaleric acid - expect True
        "OC(=O)C/C=C(\\C)/C(O)=O",                # 2-Methylglutaconic acid - expect True
        "CCCCCCCC(C)CC(O)=O",                     # 3-methylundecanoic acid - expect True
        "CCCCCCCCC(C)CCCCCC(C)CCCCCCCCCC(C)CC(O)=O",  # 3,13,19-trimethyltricosanoic acid - expect True
        "OC(C(CCCCCCCCCCCCCCCC)C)=O",             # 2-methyloctadecanoic acid - expect True
        "CC(C)CCCCCCCCCCCCCC(O)=O",               # isoheptadecanoic acid - expect True
        "CC(C)CCCCCC(O)=O",                       # 7-methyloctanoic acid - expect True
        "OC(=O)\\C=C\\C(C)C",                      # 4-Methyl-2-pentenoic acid - expect True
        "OC(=O)C(CCCC(C)C)C",                      # 2,6-dimethylheptanoic acid - expect True
        # The false positives and negatives from the previous attempt (allylglycine etc.) will now (hopefully)
        # be rejected because their overall heteroatom pattern does not meet the fatty acid criteria.
        "NC(CC=C)C(O)=O",                        # allylglycine – expect False
        "O[C@@H]([C@H](O)CCCCCCCC(O)=O)CCCCCCCC",  # threo-9,10-dihydroxystearic acid – expect False
        "OC(=O)CC(CC(O)=O)NC=O",                  # N-formylisoglutamic acid – expect False
        "CC(=CCC)C(=O)O",                         # 2-propyl-2-pentenoic acid – expect False (non–methyl branch)
    ]
    
    for smi in test_smiles_list:
        res, reason = is_methyl_branched_fatty_acid(smi)
        print(smi, res, reason)