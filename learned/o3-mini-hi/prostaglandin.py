"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Prostanoid derivatives (prostaglandins)
Definition: Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.
Heuristics used (improved):
 1. Overall carbon count roughly between 17 and 27.
 2. The molecule must contain an isolated (non-fused) cyclopentane ring consisting solely of carbons.
 3. At least one substituent branch off that ring must be “long” (≥5 carbon atoms) and contain a 
    carbonyl group motif (allowing carboxylic acid, ester, or amide groups) as a sign of derivation 
    from prostanoic acid.
If any condition is not met, the function returns False with an explanation.
"""

from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines whether the given SMILES string is likely a prostaglandin derivative.
    
    The algorithm uses the following steps:
      - Checks that the total number of carbons is in the range 17 to 27.
      - Searches for an isolated cyclopentane ring exclusively made of carbon atoms. By “isolated”
        we mean that every atom in the ring is only part of one ring (so it is not fused to other rings).
      - From atoms on that ring, we examine substituent branches (atoms attached to the ring but not
        part of it). At least one branch must have:
           * a chain length (number of carbons) of at least 5; and
           * a carbonyl type motif ([CX3](=O)[O,N]) which will catch acids, esters, or amides.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a prostaglandin derivative, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbons)
    if not (17 <= c_count <= 27):
        return False, f"Carbon count {c_count} outside prostaglandin expected range (17-27)"
    
    ring_info = mol.GetRingInfo().AtomRings()
    # Look for isolated cyclopentane rings (ring length = 5 and all atoms are carbon)
    isolated_c5_list = []
    for ring in ring_info:
        if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            # check that each atom in the ring belongs to only one ring (i.e. is not fused)
            if all(mol.GetRingInfo().NumAtomRings(idx) == 1 for idx in ring):
                isolated_c5_list.append(set(ring))
    if not isolated_c5_list:
        return False, "No isolated cyclopentane (5-carbon) ring found as prostanoid core"
    
    # Define a helper function which, given a starting atom idx (in a branch not in the core)
    # collects the substituent branch (ignoring atoms in the core ring).
    def get_branch(atom_idx, core, visited=None):
        if visited is None:
            visited = set()
        visited.add(atom_idx)
        branch_atoms = {atom_idx}
        atom = mol.GetAtomWithIdx(atom_idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in core:
                continue
            if nb_idx not in visited:
                branch_atoms.update(get_branch(nb_idx, core, visited))
        return branch_atoms

    # Build a submolecule from a set of atom indices and count how many carbons it contains.
    def count_branch_carbons(atom_ids):
        return sum(1 for idx in atom_ids if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Pattern for “carbonyl” motif that accepts acid, ester, or amide.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[O,N]")
    
    # For each isolated cyclopentane core, check the substituents attached to it.
    for core in isolated_c5_list:
        substituents = []
        # for each atom in the core, get neighbors not in the core
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx not in core:
                    substituents.append(nb_idx)
        # remove duplicates
        substituents = set(substituents)
        
        # For each substituent, obtain the entire branch (all atoms connected without crossing back into the core)
        long_branch_found = False
        carbonyl_found = False
        for sub in substituents:
            branch_atoms = get_branch(sub, core)
            branch_carbon_count = count_branch_carbons(branch_atoms)
            # Create a submol to search for the carbonyl motif
            submol = Chem.PathToSubmol(mol, list(branch_atoms))
            if submol.HasSubstructMatch(carbonyl_pattern):
                carbonyl_found = True
            # We require at least one branch with 5 or more carbon atoms
            if branch_carbon_count >= 5:
                long_branch_found = True
        if long_branch_found and carbonyl_found:
            return True, "Molecule has an isolated cyclopentane core with a long substituent that bears a carbonyl motif, and an appropriate carbon count for a prostaglandin derivative"
    
    return False, "No cyclopentane core with suitable substituents (long chain with carbonyl motif) was found"

# Uncomment the lines below for some preliminary tests:
# test_smiles = {
#    "prostaglandin A2": "CCCC[C@H](O)\\C=C\\[C@H]1C=CC(=O)[C@@H]1C\\C=C/CCCC(O)=O", 
#    "nonacosanoic acid": "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",
#    "prostaglandin F1alpha alcohol": "C(O)CCCCCC[C@@H]1[C@H]([C@@H](C[C@@H]1O)O)/C=C/[C@H](CCCCC)O",
#    "prostaglandin F2alpha dimethylamine": "CCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCCN(C)C",
# }
#
# for name, s in test_smiles.items():
#    result, reason = is_prostaglandin(s)
#    print(f"Name: {name}\nSMILES: {s}\nResult: {result}\nReason: {reason}\n")