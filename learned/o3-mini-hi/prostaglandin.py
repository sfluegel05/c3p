"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Prostanoid derivatives (prostaglandins)
Definition: Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.
Heuristics used (improved):
 1. Overall carbon count should be roughly in the range expected for prostaglandins (17 to 27 carbons).
 2. The molecule must contain an isolated (non-fused) cyclopentane ring consisting solely of carbon atoms.
 3. At least one substituent branch on that ring must be “long” (at least 5 carbons) and contain a carbonyl motif,
    as a trace of derivation from prostanoic acid.
If any condition is not met, the function returns False with an explanation.
"""

from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines whether the given SMILES string is likely a prostaglandin derivative.
    
    The algorithm uses the following steps:
      1. Attempts to parse the SMILES and, if multiple fragments exist (e.g. counter-ions present),
         selects the largest fragment.
      2. Checks that the total number of carbon atoms is in the range 17 to 27.
      3. Searches for an isolated cyclopentane ring (5-membered ring with only carbon atoms and not fused with another ring).
      4. From that cyclopentane ring, examines substituent branches to verify that at least one branch:
            - contains at least 5 carbon atoms, and
            - features a carbonyl motif ([CX3](=O)[O,N]), capturing acids, esters, or amides.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a prostaglandin derivative, False otherwise.
        str: Explanation for the decision.
    """
    # Attempt to create an RDKit molecule object
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If there are counter-ions or multiple fragments,
    # select the largest fragment (by heavy atom count).
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    except Exception:
        # In case of any error, continue with the original mol
        pass

    # Count the number of carbon atoms in the molecule
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbons)
    if not (17 <= c_count <= 27):
        return False, f"Carbon count {c_count} outside expected range for prostaglandins (17-27)"
    
    ring_info = mol.GetRingInfo().AtomRings()
    isolated_c5_list = []
    # Look for isolated cyclopentane rings: ring of length 5 in which every atom is carbon and belongs only to that ring.
    for ring in ring_info:
        if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            # Check that each atom in the ring is part of exactly 1 ring
            if all(mol.GetRingInfo().NumAtomRings(idx) == 1 for idx in ring):
                isolated_c5_list.append(set(ring))
    if not isolated_c5_list:
        return False, "No isolated cyclopentane (5-carbon) ring found as prostanoid core"
    
    # Helper function: given a starting atom not in the core, return all connected atoms (branch),
    # avoiding atoms that are in the core ring.
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

    # Count how many carbons are in a set of branch atom indices
    def count_branch_carbons(atom_ids):
        return sum(1 for idx in atom_ids if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Define the SMARTS for a carbonyl motif (captures acid, ester, amide patterns)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[O,N]")
    
    # For each isolated cyclopentane core, check substituent branches
    for core in isolated_c5_list:
        # Collect indices of atoms directly attached to the core (neighbors not in the core)
        substituent_indices = set()
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx not in core:
                    substituent_indices.add(nb_idx)
        # For each substituent, analyze its branch
        long_branch_found = False
        carbonyl_found = False
        for sub_idx in substituent_indices:
            try:
                branch_atoms = get_branch(sub_idx, core)
            except Exception as e:
                return False, f"Error processing substituent branch: {e}"
            branch_carbon_count = count_branch_carbons(branch_atoms)
            # Create a submol corresponding to the branch to search the carbonyl motif
            branch_atom_list = list(branch_atoms)
            try:
                submol = Chem.PathToSubmol(mol, branch_atom_list)
            except Exception:
                # If creation fails, skip this branch
                continue
            if submol.HasSubstructMatch(carbonyl_pattern):
                carbonyl_found = True
            if branch_carbon_count >= 5:
                long_branch_found = True
        if long_branch_found and carbonyl_found:
            return True, ("Molecule has an isolated cyclopentane core with at least one long substituent "
                          "that bears a carbonyl motif and an appropriate carbon count for a prostaglandin derivative")
    
    return False, "No cyclopentane core with suitable substituents (long branch with carbonyl motif) was found"


# Uncomment below for some preliminary tests:
# test_smiles = {
#    "prostaglandin A2": "CCCC[C@H](O)\\C=C\\[C@H]1C=CC(=O)[C@@H]1C\\C=C/CCCC(O)=O", 
#    "nonacosanoic acid": "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",
#    "prostaglandin E3": "CC\\C=C/C[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
#    "6alpha-Prostaglandin I1": "[C@@H]1([C@@H]([C@@H]2[C@H](C1)O[C@@H](C2)CCCCC(O)=O)/C=C/[C@H](CCCCC)O)O",
# }
#
# for name, s in test_smiles.items():
#    result, reason = is_prostaglandin(s)
#    print(f"Name: {name}\nSMILES: {s}\nResult: {result}\nReason: {reason}\n")