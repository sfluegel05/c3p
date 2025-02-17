"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: Diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 skeleton 
of the parent diterpene has been rearranged or modified (commonly with loss of one or more methyl groups) 
or decorated (e.g. glycosylation). Because many natural diterpenoids have additional substituents, 
heuristics based solely on the carbon count may misclassify valid examples. 

Heuristic criteria used in this function (developed iteratively based on test outcomes):
  - Only allowed elements: H, C, N, O.
  - Carbon atoms are expected to be at least 15 and no more than 50 (to catch C20-based skeletons, even if glycosylated).
  - Molecular weight should be roughly between 220 and 1200 Da.
  - A relaxed oxygen-to-carbon ratio (O/C) of no more than 0.6 is allowed (some glycosides are more oxygenated).
  - Nitrogen is not common in diterpenoids so molecules with more than one N are rejected.
  - Diterpenoids usually arise from four isoprene units and are not simple linear fatty acid chains;
    thus if an acyclic molecule is completely linear (i.e. the longest contiguous carbon‐chain equals total carbon count),
    it is likely a fatty acid derivative rather than a diterpenoid.
  - A high fraction of aromatic carbons is not typical.
If these criteria are met, the molecule is considered likely to be a diterpenoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest chain consisting solely of carbon atoms.
    Builds an undirected graph of carbon atoms and performs a depth-first search from each carbon.
    """
    # Get indices for carbon atoms only.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return 0
    # Build an adjacency list of carbon neighbors.
    neighbors = {idx: [] for idx in carbon_idxs}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
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
    Determines if a molecule is likely a diterpenoid based on its SMILES string using improved heuristics.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is likely a diterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Allowed Elements.
    allowed_atomic_nums = {1, 6, 7, 8}  # H, C, N, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} which is not typical for diterpenoids."
            
    # Count different atoms.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    nO = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    nN = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Criterion 2: Carbon count (relaxed upper bound to allow glycosides etc).
    if nC < 15:
        return False, f"Too few carbon atoms ({nC}) for a diterpenoid skeleton."
    if nC > 50:
        return False, f"Too many carbon atoms ({nC}); exceeds expected diterpenoid range."
    
    # Criterion 3: Molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 220 or mw > 1200:
        return False, f"Molecular weight {mw:.1f} Da is outside the typical range (220–1200 Da) for diterpenoids."
        
    # Criterion 4: Oxygen-to-carbon ratio.
    o_c_ratio = nO / nC
    if o_c_ratio > 0.6:
        return False, f"O/C ratio is too high ({o_c_ratio:.2f}); diterpenoids are usually modestly oxygenated."
        
    # Criterion 5: Nitrogen count.
    if nN > 1:
        return False, f"Too many nitrogen atoms ({nN}) for a typical diterpenoid."
        
    # Criterion 6: Aromaticity.
    aromatic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    aromatic_ratio = aromatic_carbons / nC
    if aromatic_ratio > 0.5:
        return False, f"High aromaticity detected (aromatic C fraction: {aromatic_ratio:.2f}); diterpenoids are usually aliphatic or partly aliphatic."
        
    # Criterion 7: Ring system or carbon chain branching.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings > 0:
        ring_msg = f"{n_rings} ring(s) detected."
    else:
        # For acyclic molecules, check if the longest carbon-only chain runs through all carbons
        lc = longest_carbon_chain(mol)
        if lc == nC:
            return False, f"Acyclic and appears linear (longest carbon chain = {lc} equals total C count = {nC}); likely a fatty acid derivative."
        if lc < 10:
            return False, f"Acyclic and longest carbon chain is only {lc} atoms long; not consistent with a diterpenoid skeleton."
        ring_msg = f"Acyclic with longest carbon chain of {lc} atoms."
    
    return True, (f"Likely diterpenoid: {nC} carbons, {mw:.1f} Da, {ring_msg} "
                  f"O/C ratio = {o_c_ratio:.2f}, aromatic C fraction = {aromatic_ratio:.2f}, N count = {nN}")

# Example usage:
if __name__ == "__main__":
    # Test with gibberellin A34 as an example diterpenoid.
    test_smiles = "[H][C@]12CC[C@]3([H])[C@](CC1=C)(C2)[C@@H](C(O)=O)[C@]1([H])[C@@]2(C)[C@@H](O)[C@@H](O)C[C@@]31OC2=O"
    result, explanation = is_diterpenoid(test_smiles)
    print(result, explanation)