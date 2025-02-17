"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: Sphingoid compounds
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.
Improved strategy:
  - Instead of checking whole disconnected fragments, we search the entire molecule for an acyclic sphingoid core
    (using several SMARTS patterns).
  - For every match we extract a neighborhood (all atoms within a given bond distance from the core) to ensure
    that the long aliphatic chain and the amino group (i.e. nitrogen) are part of the same connected region.
  - We then require that this “sphingoid fragment” has:
       * A long aliphatic chain (at least 8 contiguous non‐ring carbons),
       * At least one nitrogen atom,
       * And that the core match atoms are all acyclic.
  - If any neighborhood meets these criteria, the molecule is classified as sphingoid.
  
This approach helps to avoid false positives where the two features are separated in different parts of the molecule.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid compound based on its SMILES string.
    Sphingoid compounds include sphinganine (and its homologs/stereoisomers) and the hydroxy/unsaturated variants.
    
    For classification we:
      - Search for an acyclic sphingoid core using several SMARTS patterns.
      - For each core match, build a connected “neighborhood” (within a specified bond radius) and 
        require that it contains:
           * A long aliphatic chain (at least 8 contiguous, non‐ring carbons, with single or double bonds),
           * At least one nitrogen atom.
      - Basic fragment properties (e.g. approximate molecular weight, carbon and nitrogen counts) are used
        to aid in explanation.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): A tuple containing the classification (True if sphingoid is detected, False otherwise)
                     and an explanation message.
    """
    # Parse the input molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a relaxed SMARTS pattern for a long aliphatic chain:
    # Matches 8 contiguous non‐ring carbon atoms (allowing either sp3 or sp2 atoms).
    chain_smarts = "[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    if chain_pattern is None:
        return False, "Error in chain SMARTS pattern"
    
    # Define several sphingoid core SMARTS patterns.
    # These patterns capture variants like a fully hydroxylated amino-diol core, a carbonyl variant etc.
    core_smarts_list = [
        "C(O)C(N)CO",    # fully hydroxylated amino-diol core
        "C(=O)C(N)CO",   # carbonyl variant (dehydro form)
        "C(=O)CN",       # deoxy-carbonyl variant
        "C(O)CN"         # deoxy-hydroxy variant
    ]
    core_patterns = []
    for s in core_smarts_list:
        pat = Chem.MolFromSmarts(s)
        if pat is not None:
            core_patterns.append(pat)

    # Helper function: return True if all atoms in a given match (by indices) are acyclic.
    def acyclic_match(mol_obj, match_indices):
        return all(not mol_obj.GetAtomWithIdx(idx).IsInRing() for idx in match_indices)
    
    # Helper function: given a set of seed atom indices, return all atoms within 'radius' bonds.
    def get_neighborhood_atoms(mol_obj, seed_indices, radius):
        visited = set(seed_indices)
        frontier = set(seed_indices)
        for _ in range(radius):
            next_frontier = set()
            for idx in frontier:
                atom = mol_obj.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx not in visited:
                        visited.add(n_idx)
                        next_frontier.add(n_idx)
            frontier = next_frontier
        return visited

    # Loop over each sphingoid core pattern.
    for pat in core_patterns:
        matches = mol.GetSubstructMatches(pat)
        for match in matches:
            # Only consider the match if all atoms in it are acyclic (not in any ring)
            if not acyclic_match(mol, match):
                continue
            # Build a neighborhood around the core match.
            # We use a radius of 4 bonds (this parameter can be tuned).
            neighborhood = get_neighborhood_atoms(mol, list(match), radius=4)
            # Create a submolecule from the neighborhood atoms.
            submol = Chem.PathToSubmol(mol, list(neighborhood))
            # Check: The submol must contain a long aliphatic chain.
            if not submol.HasSubstructMatch(chain_pattern):
                continue
            # Also require that the submol has at least one nitrogen atom (amino group)
            n_count = sum(1 for atom in submol.GetAtoms() if atom.GetAtomicNum() == 7)
            if n_count < 1:
                continue

            # (Optional) Check additional properties. Here we compute weight and carbon count.
            frag_wt = rdMolDescriptors.CalcExactMolWt(submol)
            c_count = sum(1 for atom in submol.GetAtoms() if atom.GetAtomicNum() == 6)
            # Accept a broad range here; sphingoid fragments usually lie between 200 and 1000 Da 
            # and contain a moderate-to-high carbon count.
            if frag_wt < 200 or frag_wt > 1000 or c_count < 16:
                continue

            # If all tests are passed, we consider the molecule a sphingoid compound.
            reason = ("Fragment with weight {:.1f} Da, {} carbons, and {} nitrogen(s) contains an acyclic sphingoid core"
                      " and a long aliphatic chain within a {}-bond neighborhood.".format(frag_wt, c_count, n_count, 4))
            return True, reason

    return False, "No connected substructure combining an acyclic sphingoid core, long aliphatic chain, and amino group was found."

# For testing purposes, run a few examples.
if __name__ == "__main__":
    test_smiles = [
        # True positives.
        "CCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C",  # Example glycosylated sphingoid
        "CCCCCCCCCCCC\\C=C\\C(=O)[C@@H](N)CO",  # 3-dehydrosphingosine
        "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C",  # N-heptadecanoyl-... sphingoid
        "CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)N(C)C",  # N,N-dimethylsphingosine
        # A false positive candidate that should not classify as sphingoid.
        "CCCC\\C=C/C\\C=C/CCCCCCCCCC(=O)NCCO"
    ]
    for smi in test_smiles:
        result, reason = is_sphingoid(smi)
        print("SMILES:", smi)
        print("Is sphingoid?:", result)
        print("Reason:", reason)
        print()