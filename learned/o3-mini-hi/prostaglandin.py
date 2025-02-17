"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Prostanoid derivatives (prostaglandins)
Definition: Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.
Heuristics (improved):
 1. Total carbon atom count should be roughly in the range expected for prostaglandins (17-30).
 2. The molecule must contain at least one 5-membered ring composed exclusively of carbon atoms.
 3. The molecule must have a carbonyl as part of an acid/ester motif (SMARTS: C(=O)O).
 4. At least one substituent branch off the cyclopentane ring must be “long” (≥5 carbons) and contain a carbonyl (SMARTS: [CX3](=O)).
If any condition is not met, the function returns False with an explanation.
"""

from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines whether the given SMILES string is likely a prostaglandin derivative.
    
    The algorithm uses the following steps:
      1. Parse the SMILES and if multiple fragments exist, select the largest.
      2. Count the number of carbon atoms; if the count is not in the range 17-30, assume it is not prostaglandin.
      3. Look for at least one 5-membered ring that is composed entirely of carbon atoms.
      4. Verify that the molecule contains an acid/ester carbonyl motif (C(=O)O).
      5. For at least one candidate cyclopentane ring, examine all substituent branches (neighbors not in the ring)
         and recursively collect the branch atoms. If at least one branch contains at least 5 carbons and a carbonyl (using SMARTS [CX3](=O)), 
         classify the molecule as a prostaglandin derivative.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a prostaglandin derivative, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If multiple fragments exist, select the largest by heavy atom count.
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    except Exception:
        pass
    
    # Count total carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbons)
    if not (17 <= c_count <= 30):
        return False, f"Carbon count {c_count} outside expected range for prostaglandins (17-30)"
    
    # Check if the molecule has an acid/ester carbonyl motif.
    # (This SMARTS should capture carboxylic acids and simple esters.)
    acid_ester_smarts = "C(=O)O"
    acid_ester_pattern = Chem.MolFromSmarts(acid_ester_smarts)
    if not mol.HasSubstructMatch(acid_ester_pattern):
        return False, "No acid/ester carbonyl motif (C(=O)O) found"
    
    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    cyclopentane_rings = []
    for ring in ring_info:
        if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            cyclopentane_rings.append(set(ring))
    if not cyclopentane_rings:
        return False, "No 5-membered cyclopentane ring (all carbon atoms) found as prostanoid core"
    
    # Define SMARTS for a standalone carbonyl (e.g. part of a ketone, acid, or ester).
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)")
    
    # Recursive helper function: given a starting atom (by index) not in the core, get the full connected branch.
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
    
    # For each candidate cyclopentane ring, check its substituent branches.
    for core in cyclopentane_rings:
        substituent_indices = set()
        # For every atom in the ring, collect neighbors not in the ring.
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx not in core:
                    substituent_indices.add(nb_idx)
        long_branch_with_carbonyl_found = False
        for sub_idx in substituent_indices:
            try:
                branch_atom_ids = get_branch(sub_idx, core)
            except Exception as e:
                return False, f"Error processing substituent branch: {e}"
            # Count carbons in the branch.
            branch_carbons = sum(1 for idx in branch_atom_ids if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # Create a submolecule from branch atoms for substructure search.
            try:
                branch_mol = Chem.PathToSubmol(mol, list(branch_atom_ids))
            except Exception:
                continue
            if branch_carbons >= 5 and branch_mol.HasSubstructMatch(carbonyl_pattern):
                long_branch_with_carbonyl_found = True
                break
        if long_branch_with_carbonyl_found:
            return True, ("Molecule has a 5-membered carbocycle combined with a long substituent branch "
                          "bearing a carbonyl motif (as expected for a prostanoid derivative) and a suitable carbon count.")
    
    return False, "No cyclopentane core with a long substituent branch containing a carbonyl motif was found"

# Example usage:
if __name__ == "__main__":
    test_examples = {
        "prostaglandin A2": "CCCC[C@H](O)\\C=C\\[C@H]1C=CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
        "nonacosanoic acid": "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # Fatty acid, should be false.
        "prostaglandin E3": "CC\\C=C/C[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
    }
    for name, s in test_examples.items():
        result, reason = is_prostaglandin(s)
        print(f"Name: {name}\nSMILES: {s}\nResult: {result}\nReason: {reason}\n")