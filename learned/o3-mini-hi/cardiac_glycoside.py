"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside
Definition: Steroid lactones containing sugar residues that act on the contractile force 
            of the cardiac muscles.
This program checks for the presence of:
  1) A steroid nucleus (a fused tetracyclic ring system: roughly 3 six‐membered rings and 1 five‐membered ring)
  2) A lactone ring (approximated as a butenolide / furan-2-one substructure)
  3) At least one sugar residue (detected via a pyranose-like ring pattern, with a looser fallback)
Note: The SMARTS definitions and fused ring detection are approximate and may not capture every nuance.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Cardiac glycosides are steroid lactones with sugar moieties.
    This function checks for:
      - A steroid nucleus (fused tetracyclic system, using ring-fusion heuristics)
      - A lactone group (butenolide/furan-2-one)
      - At least one sugar residue (pyranose-like ring)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a cardiac glycoside, False otherwise.
        str: Reason for the classification outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Check for a steroid nucleus (fused tetracyclic system)
    # We use the heuristic that a steroid nucleus should contain at least 4 fused rings -
    # with rings of size 5 or 6. Fused rings are defined as rings sharing at least 2 atoms.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    if len(candidate_rings) < 4:
        return False, "Less than 4 rings of size 5 or 6 found; no steroid nucleus detected"
    
    # Build a simple graph: nodes are ring indices; add an edge if two rings share at least 2 atoms.
    adj = {i: set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    # Find connected components (via DFS)
    visited = set()
    def dfs(node, comp):
        comp.add(node)
        visited.add(node)
        for nbr in adj[node]:
            if nbr not in visited:
                dfs(nbr, comp)
        return comp
    
    found_steroid = False
    for i in range(len(candidate_rings)):
        if i not in visited:
            comp = dfs(i, set())
            if len(comp) >= 4:
                found_steroid = True
                break
    if not found_steroid:
        return False, "No fused tetracyclic steroid nucleus detected"

    # (2) Check for a lactone group: use a SMARTS for butenolide (furan-2-one)
    lactone_smarts = "O=C1C=CCO1"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group (butenolide/furan-2-one) found"
    
    # (3) Check for at least one sugar residue.
    # Primary try: a stereochemically defined pyranose ring
    sugar_smarts = "[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if mol.HasSubstructMatch(sugar_pattern):
        sugar_found = True
    else:
        # Fallback: a looser pattern for a six-membered ring with at least two oxygens.
        sugar_smarts_loose = "O1CC(O)CC(O)C1"
        sugar_pattern_loose = Chem.MolFromSmarts(sugar_smarts_loose)
        sugar_found = mol.HasSubstructMatch(sugar_pattern_loose)
    if not sugar_found:
        return False, "No sugar residue (pyranose-like ring) found"

    # Optionally, check that the molecule's molecular weight is not too low
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for typical cardiac glycosides"

    return True, "Contains a fused steroid nucleus, a lactone ring, and at least one sugar residue"