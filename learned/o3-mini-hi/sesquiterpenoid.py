"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: Sesquiterpenoid
Definition:
  Any terpenoid derived from a sesquiterpene (C15 skeleton) even if the skeleton
  has been rearranged or modified (e.g., by removal of one or more skeletal atoms).
  Our heuristic now checks:
    - Allowed elements: C, H, O, and S.
    - Carbon count between 12 and 21.
    - Molecular weight between 150 and 500 Da.
    - Aromatic atoms are permitted if they do not form a high fraction (>50%) of total carbons.
    - For acyclic molecules, if the structure is highly flexible (rotatable bonds >~6)
      or has a nearly unbranched carbon chain together with a carboxylic acid group,
      then it is likely a fatty acid rather than a sesquiterpenoid.
Note: This heuristic is approximate and may not cover every edge case.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def _longest_carbon_chain(mol):
    """
    For an acyclic molecule, computes the length (number of carbon atoms) of the longest
    continuous chain (considering only atoms with atomic number 6).
    This is done by a simple DFS on the carbon-only subgraph.
    """
    # Build a graph: node indices for carbons only.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Create a dictionary of carbon neighbors (only consider carbons)
    neighbors = {}
    for idx in carbon_idxs:
        nbrs = []
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                nbrs.append(nbr.GetIdx())
        neighbors[idx] = nbrs

    longest = 0

    def dfs(current, visited):
        max_length = 1
        for nbr in neighbors.get(current, []):
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    for idx in carbon_idxs:
        chain_length = dfs(idx, {idx})
        if chain_length > longest:
            longest = chain_length
    return longest

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    
    Our improved heuristic requires:
      (1) Allowed elements: only hydrogen, carbon, oxygen, and sulfur.
      (2) Carbon count between 12 and 21.
      (3) Molecular weight between 150 and 500 Da.
      (4) Aromatic atoms are permitted provided that the fraction of aromatic carbons 
          is not more than 50%.
      (5) For acyclic molecules (no rings), if the molecule is very flexible (many rotatable bonds)
          or if it is an unbranched (linear) structure bearing a carboxylic acid group,
          then it is likely a fatty acid rather than a sesquiterpenoid.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a sesquiterpenoid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check allowed elements:
    # Allow only H, C, O, S (atomic numbers 1,6,8,16).
    allowed_atomic_nums = {1, 6, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()}, not allowed for a typical sesquiterpenoid"
    
    # 2. Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (12 <= c_count <= 21):
        return False, f"Molecule contains {c_count} carbons, not consistent with a typical sesquiterpene-derived structure (expected roughly 12–21 carbons)"
    
    # 3. Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 500):
        return False, f"Molecular weight {mol_wt:.1f} Da not in expected range (150–500 Da) for sesquiterpenoids"
    
    # 4. Check aromaticity:
    # Instead of rejecting any aromatic atoms, allow if the fraction of aromatic carbons is modest.
    aromatic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    if c_count > 0 and (aromatic_carbons / c_count) > 0.5:
        return False, "Too many aromatic carbons (more than 50%), which is not typical for a sesquiterpenoid structure"
    
    # 5. Extra checks for acyclic molecules.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if num_rings == 0:
        # (a) If the molecule is acyclic and very flexible, it might be a fatty acid.
        if n_rotatable > 6:
            return False, f"Acyclic structure with high flexibility (rotatable bonds = {n_rotatable}), likely a fatty acid rather than a sesquiterpenoid"
        
        # (b) If an acyclic molecule bears a carboxylic acid group and shows a nearly linear,
        #     unbranched carbon chain, then it is likely a fatty acid.
        acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
        if mol.HasSubstructMatch(acid_smarts):
            longest_chain = _longest_carbon_chain(mol)
            if longest_chain >= (c_count - 2):
                return False, f"Acyclic fatty acid: longest carbon chain length ({longest_chain}) is almost the entire carbon count ({c_count}) in a molecule bearing a carboxylic acid"
    
    return True, f"Molecule has {c_count} carbons, MW {mol_wt:.1f} Da, and passes allowed element and structural checks, consistent with a sesquiterpene-derived structure"

# Example usage (for testing purposes)
if __name__ == '__main__':
    # Test with a few examples:
    tests = {
        "indicanone": "CC(=C)[C@@H]1CCC(CO)=C2CC(=O)C(C)=C2C1",
        "juvenile hormone I acid": "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",
        "(-)-GR24": "O1C(/C(/[C@H]2[C@@H]1C3=C(C2)C=CC=C3)=C/O[C@H]4OC(C(=C4)C)=O)=O",
        "3,7,11-trimethyl-dodecanoic acid": "OC(=O)CC(CCCC(CCCC(C)C)C)",
        "geranyl acetone (false positive expected)": "CC(=O)CC\\C=C(/C)CCC=C(C)C",
    }
    for name, smi in tests.items():
        result, reason = is_sesquiterpenoid(smi)
        print(f"{name}: {result} -- {reason}")