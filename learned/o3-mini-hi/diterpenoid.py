"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: Diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 skeleton 
of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms 
(generally methyl groups).

Heuristic criteria used in this function:
  - The elements should be mostly C, H, O, and optionally N.
  - The number of carbon atoms is expected to be roughly between 15 and 33.
  - The molecular weight is expected to be in the range 220–800 Da.
  - The oxygen-to-carbon ratio (O/C) should be low (≤0.5) because diterpenoids are not highly oxygenated.
  - Many diterpenoids are cyclic; however, for acyclic molecules we additionally check that the carbon 
    skeleton is not completely linear. If an acyclic molecule’s longest strictly carbon-only chain equals 
    the total carbon count it likely comes from a fatty acid.
If these criteria are met, the molecule is considered likely to be a diterpenoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest chain that consists only of carbon atoms.
    This is done by building an undirected graph over carbons and doing a DFS search on each node.
    For small molecules the computational cost is acceptable.
    """
    # Get indices for carbon atoms
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return 0
    # Build an adjacency list for carbon atoms only.
    neighbors = {idx: [] for idx in carbon_idxs}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # only consider carbon neighbors
                neighbors[idx].append(nbr.GetIdx())
    
    best = 0
    def dfs(current, visited):
        nonlocal best
        best = max(best, len(visited))
        for nbr in neighbors[current]:
            if nbr not in visited:
                dfs(nbr, visited | {nbr})
    
    for start in carbon_idxs:
        dfs(start, {start})
    return best

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string using heuristic checks.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if molecule is likely a diterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Criterion 1: Allowed Elements.
    # We allow hydrogen (1), carbon (6), oxygen (8) and nitrogen (7) because some diterpenoids include N.
    allowed_atomic_nums = {1, 6, 7, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} which is not typical for diterpenoids."
    
    # Criterion 2: Carbon count.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 15:
        return False, f"Too few carbon atoms ({nC}) for a diterpenoid skeleton."
    if nC > 33:
        return False, f"Too many carbon atoms ({nC}); exceeds expected diterpenoid range."
    
    # Criterion 3: Molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 220 or mw > 800:
        return False, f"Molecular weight {mw:.1f} Da is outside the typical range (220-800 Da) for diterpenoids."
    
    # Criterion 4: Oxygen-to-Carbon ratio.
    nO = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    o_c_ratio = nO / nC
    if o_c_ratio > 0.5:
        return False, f"O/C ratio is too high ({o_c_ratio:.2f}); diterpenoids are typically less oxygenated."
    
    # Criterion 5: Ring system or carbon chain branching.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings > 0:
        ring_msg = f"{n_rings} ring(s) detected."
    else:
        # For acyclic molecules we check the longest carbon chain.
        lc = longest_carbon_chain(mol)
        # If the longest chain equals the total carbon count, the structure is very linear (like a fatty acid).
        if lc == nC:
            return False, f"Acyclic and appears linear ({lc} C atoms in longest chain equals total C count); likely a fatty acid."
        if lc < 10:
            return False, f"No cyclic ring system and longest carbon chain is only {lc} atoms long; diterpenoids typically have a ~C20 backbone."
        ring_msg = f"Acyclic with longest carbon chain of {lc} atoms."
    
    # All criteria are met.
    return True, f"Likely diterpenoid: {nC} carbons, {mw:.1f} Da, {ring_msg} O/C ratio = {o_c_ratio:.2f}"

# Example usage (can be removed in production):
if __name__ == "__main__":
    # Try one provided example: gibberellin A34.
    test_smiles = "[H][C@]12CC[C@]3([H])[C@](CC1=C)(C2)[C@@H](C(O)=O)[C@]1([H])[C@@]2(C)[C@@H](O)[C@@H](O)C[C@@]31OC2=O"
    result, reason = is_diterpenoid(test_smiles)
    print(result, reason)