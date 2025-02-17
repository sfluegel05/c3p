"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside
Definition: Steroid lactones containing sugar residues that act on the contractile force
            of the cardiac muscles.
This program checks for the presence of:
  1) A steroid nucleus (a fused tetracyclic ring system – roughly 3 six‐membered 
     rings and 1 five‐membered ring fused together).
  2) A lactone ring (approximated as a butenolide/furan-2-one substructure).
  3) At least one sugar residue (detected by scanning for a pyranose/furanose‐like ring, 
     i.e. a 6-membered ring with at least 2 oxygens or a 5-membered ring with at least 1 oxygen,
     that is not “buried” within the steroid core).
Note: The SMARTS definitions and fused‐ring detection are approximations.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Cardiac glycosides are steroid lactones with sugar moieties.
    This function checks for:
      - A steroid nucleus: detected by looking for a fused system of at least 4 rings 
        (of size 5 or 6) that are bonded via at least 2 overlapping atoms.
      - A lactone ring: matched via a butenolide/furan-2-one SMARTS pattern.
      - At least one sugar residue: found by scanning all rings for a pyranose (6-membered ring with >=2 O)
        or furanose (5-membered ring with >=1 O) that is attached (sharing at most 1 atom) to the steroid nucleus.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a cardiac glycoside, False otherwise.
        str: Reason for the classification outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Detect steroid nucleus:
    # First, get all rings of size 5 or 6.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    if len(candidate_rings) < 4:
        return False, "Less than 4 rings of size 5 or 6 found; no steroid nucleus detected"
    
    # Build a simple graph where each node corresponds to a candidate ring.
    # Two nodes are connected if they share at least 2 atoms.
    adj = {i: set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
                
    # Find connected components (each is a set of indices referring to candidate_rings)
    visited = set()
    def dfs(node, comp):
        comp.add(node)
        visited.add(node)
        for nbr in adj[node]:
            if nbr not in visited:
                dfs(nbr, comp)
        return comp
    
    steroid_found = False
    steroid_atoms = set()
    for i in range(len(candidate_rings)):
        if i not in visited:
            comp = dfs(i, set())
            if len(comp) >= 4:
                steroid_found = True
                # The steroid nucleus is (approximately) the union of all atoms in these rings.
                for idx in comp:
                    steroid_atoms.update(candidate_rings[idx])
                break
    if not steroid_found:
        return False, "No fused tetracyclic steroid nucleus detected"
    
    # (2) Check for a lactone group. We use a SMARTS for a butenolide (furan-2-one).
    lactone_smarts = "O=C1C=CCO1"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group (butenolide/furan-2-one) found"
    
    # (3) Check for at least one sugar residue.
    # We scan all rings from the molecule (using ring info) and look for rings that look like sugars.
    sugar_found = False
    for ring in ring_info:
        ring_set = set(ring)
        # Skip rings that are essentially part of the steroid nucleus. 
        # We allow if the overlap with steroid nucleus is <= 1.
        overlap = len(ring_set.intersection(steroid_atoms))
        # Consider rings of size 6 (pyranose-like) and size 5 (furanose-like).
        if len(ring) == 6:
            # Count number of oxygen atoms in the ring.
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count >= 2 and overlap <= 1:
                sugar_found = True
                break
        elif len(ring) == 5:
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count >= 1 and overlap <= 1:
                sugar_found = True
                break
    # Fallback: if not found via ring analysis, try a loose SMARTS match.
    if not sugar_found:
        sugar_smarts_loose = "O1CC(O)CC(O)C1"
        sugar_pattern_loose = Chem.MolFromSmarts(sugar_smarts_loose)
        if mol.HasSubstructMatch(sugar_pattern_loose):
            sugar_found = True
            
    if not sugar_found:
        return False, "No sugar residue (pyranose/furanose-like ring) found"
    
    # Optionally, add a molecular weight filter.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for typical cardiac glycosides"
    
    return True, "Contains a fused steroid nucleus, a lactone ring, and at least one sugar residue"

# (Optional) Example usage:
if __name__ == '__main__':
    test_smiles = [
        # Example: gitaloxin (should be cardiac glycoside; might have borderline sugar detection)
        "[H]C(=O)O[C@H]1C[C@]2(O)[C@]3([H])CC[C@]4([H])C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]2(C)[C@@]1([H])C1=CC(=O)OC1)"
    ]
    for smi in test_smiles:
        res, reason = is_cardiac_glycoside(smi)
        print("Result:", res, "| Reason:", reason)