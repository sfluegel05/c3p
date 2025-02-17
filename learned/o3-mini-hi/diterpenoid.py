"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: Diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 skeleton 
of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms 
(generally methyl groups).

Heuristic criteria used in this function:
  - Only allowed elements: H, C, N, O.
  - Carbon atoms should be roughly between 15 and 30.
  - Molecular weight in the range 220–800 Da.
  - Oxygen-to-carbon ratio (O/C) should be low (≤0.5) since diterpenoids are not extensively oxygenated.
  - Diterpenoids are usually built around an isoprene (C5) backbone so too high aromaticity is unusual.
  - For acyclic molecules, if the longest connected carbon-only chain equals the total carbon count,
    the molecule is likely a fatty acid derivative.
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
    Determines if a molecule is a diterpenoid based on its SMILES string using improved heuristic checks.
    
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
    
    # Criterion 2: Carbon count - expect roughly from 15 to 30 carbons.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 15:
        return False, f"Too few carbon atoms ({nC}) for a diterpenoid skeleton."
    if nC > 30:
        return False, f"Too many carbon atoms ({nC}); exceeds expected diterpenoid range."
    
    # Criterion 3: Molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 220 or mw > 800:
        return False, f"Molecular weight {mw:.1f} Da is outside the typical range (220–800 Da) for diterpenoids."
    
    # Criterion 4: Oxygen-to-carbon ratio.
    nO = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    o_c_ratio = nO / nC
    if o_c_ratio > 0.5:
        return False, f"O/C ratio is too high ({o_c_ratio:.2f}); diterpenoids are typically less oxygenated."
    
    # Criterion 5: Check for over-aromatization.
    aromatic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    aromatic_ratio = aromatic_carbons / nC
    if aromatic_ratio > 0.5:
        return False, f"High aromaticity detected (aromatic C fraction: {aromatic_ratio:.2f}); diterpenoids are usually aliphatic or partly aliphatic."
    
    # Criterion 6: Ring system or carbon chain branching.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings > 0:
        ring_msg = f"{n_rings} ring(s) detected."
    else:
        # For acyclic molecules, check that the structure is not completely linear.
        lc = longest_carbon_chain(mol)
        if lc == nC:
            return False, f"Acyclic and appears linear (longest carbon chain = {lc} equals total C count = {nC}); likely a fatty acid derivative."
        if lc < 10:
            return False, f"Acyclic and longest carbon chain is only {lc} atoms long; typical diterpenoids have a more branched C20 backbone."
        ring_msg = f"Acyclic with longest carbon chain of {lc} atoms."
    
    # Optionally, one could check for isoprene-like substructures.
    # isoprene_unit = Chem.MolFromSmarts("C(C)=CC")
    # has_isoprene = mol.HasSubstructMatch(isoprene_unit)
    # This check is not mandatory because rearrangements may obscure the pattern.
    
    return True, f"Likely diterpenoid: {nC} carbons, {mw:.1f} Da, {ring_msg} O/C ratio = {o_c_ratio:.2f}, aromatic C fraction = {aromatic_ratio:.2f}"

# Example usage (can be removed in production):
if __name__ == "__main__":
    # Test with gibberellin A34.
    test_smiles = "[H][C@]12CC[C@]3([H])[C@](CC1=C)(C2)[C@@H](C(O)=O)[C@]1([H])[C@@]2(C)[C@@H](O)[C@@H](O)C[C@@]31OC2=O"
    result, reason = is_diterpenoid(test_smiles)
    print(result, reason)