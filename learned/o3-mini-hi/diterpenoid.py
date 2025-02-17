"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: Diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 skeleton 
of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms 
(generally methyl groups).

Heuristic criteria used in this function:
  - The elements should be mostly C, H, and O.
  - The number of carbon atoms is expected to be roughly between 15 and 33.
  - The molecular weight is expected to be in the range 220–800 Da.
  - The oxygen-to-carbon ratio (O/C) is expected to be low (≤0.5) because diterpenoids are not as oxygen‐rich as e.g. carbohydrates.
  - Many diterpenoids are cyclic; however, in order not to miss acyclic diterpenoids we additionally check
    that in acyclic molecules the longest connected chain of carbon atoms is long (≥10 atoms) to capture the four isoprene units.
If these criteria are met, the molecule is considered likely to be a diterpenoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest chain consisting only of carbon atoms.
    This is done by building an undirected graph of carbon atoms and doing a DFS search on each.
    (For small molecules the exponential search is acceptable.)
    """
    # Get list of indices of carbon atoms
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return 0
    # Build adjacency list only for carbon atoms
    neighbors = {idx: [] for idx in carbon_idxs}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # only consider carbon neighbors
                neighbors[idx].append(nbr.GetIdx())
    
    # DFS to find longest simple path in the carbon subgraph.
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
    Determines if a molecule is a diterpenoid based on its SMILES string using an improved set of heuristic checks.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a diterpenoid, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Criterion 1: Allowed Elements. We require that every atom be one of C, H, or O.
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} which is not typical for diterpenoids."
    
    # Criterion 2: Carbon count.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 15:
        return False, f"Too few carbon atoms ({nC}) for a diterpenoid skeleton."
    if nC > 33:
        return False, f"Too many carbon atoms ({nC}); exceeds expected diterpenoid range."
    
    # Criterion 3: Molecular weight check.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 220 or mw > 800:
        return False, f"Molecular weight {mw:.1f} Da is outside the typical range (220-800 Da) for diterpenoids."
    
    # Criterion 4: Oxygen-to-Carbon ratio.
    nO = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    o_c_ratio = nO / nC
    if o_c_ratio > 0.5:
        return False, f"O/C ratio is too high ({o_c_ratio:.2f}); diterpenoids are typically less oxygenated."
    
    # Criterion 5: Ring system or chain length.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings > 0:
        # (If rings exist we do not insist on non-aromatic rings only – many diterpenoids have fused rings.)
        ring_msg = f"{n_rings} ring(s) detected."
    else:
        # For acyclic molecules, compute the longest carbon chain as a proxy for the terpene backbone.
        lc = longest_carbon_chain(mol)
        if lc < 10:
            return False, f"No cyclic ring system and longest carbon chain is only {lc} atoms long; diterpenoids typically have a C20 (or close) backbone."
        ring_msg = f"Acyclic with longest carbon chain of {lc} atoms."
    
    # If all criteria are passed, the molecule is considered likely a diterpenoid.
    return True, (f"Likely diterpenoid: {nC} carbons, {mw:.1f} Da, {ring_msg} O/C ratio = {o_c_ratio:.2f}")

# Example use (these lines can be removed or commented out in production):
if __name__ == "__main__":
    # Try one of the provided examples: gibberellin A34.
    test_smiles = "[H][C@]12CC[C@]3([H])[C@](CC1=C)(C2)[C@@H](C(O)=O)[C@]1([H])[C@@]2(C)[C@@H](O)[C@@H](O)C[C@@]31OC2=O"
    result, reason = is_diterpenoid(test_smiles)
    print(result, reason)