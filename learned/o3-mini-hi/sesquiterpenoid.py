"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: Sesquiterpenoid
Definition:
  Any terpenoid derived from a sesquiterpene (C15 skeleton) even if the skeleton
  has been rearranged or modified (e.g., by removal of one or more skeletal atoms).
  
Heuristic improvements:
  - Allowed elements: only hydrogen, carbon, oxygen, and sulfur.
  - Carbon count is now required to be between 14 and 21 (instead of 12–21).
  - Molecular weight must be between 150 and 500 Da.
  - The fraction of aromatic carbons must be ≤50%.
  - For acyclic molecules: if the structure is very flexible (rotatable bonds >6) or
    bears a carboxylic acid group in an almost linear carbon chain, then it is likely an 
    acyclic fatty acid rather than a sesquiterpenoid.
  - Additionally, we check that the molecule shows at least three alkyl methyl groups 
    (CH3 groups attached to another carbon rather than to oxygen or sulfur).
Note: This heuristic is approximate and may still miss edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def _longest_carbon_chain(mol):
    """
    For an acyclic molecule, compute the length (number of carbon atoms) of the longest
    continuous chain (considering only atoms with atomic number 6).
    This uses a simple depth-first search on the subgraph of carbons.
    """
    # Only carbon atoms
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Build neighbors for carbons only.
    neighbors = {}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        nbrs = []
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
    Determines whether a molecule is a sesquiterpenoid (or derived from a sesquiterpene) using heuristic rules.
    
    The following checks are used:
      1. Allowed elements: only H, C, O, and S.
      2. Carbon count must be roughly 14–21.
      3. Molecular weight must lie between 150 and 500 Da.
      4. The fraction of aromatic carbons (if any) must be ≤50%.
      5. If the molecule is acyclic (no rings) and either very flexible (>6 rotatable bonds)
         or has an almost unbranched carbon chain with a carboxylic acid group, then it is considered
         likely a fatty acid rather than a sesquiterpenoid.
      6. Additionally, the molecule should have at least three “alkyl” methyl groups (CH₃ groups bonded to carbon
         rather than to an electronegative atom) – a typical feature of terpene-derived structures.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a sesquiterpenoid.
      str: A message explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Allowed elements: only H, C, O, and S.
    allowed_atomic_nums = {1, 6, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()}, not allowed for a typical sesquiterpenoid"
    
    # 2. Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (14 <= c_count <= 21):
        return False, f"Molecule contains {c_count} carbons, not consistent with a typical sesquiterpene-derived structure (expected roughly 14–21 carbons)"
    
    # 3. Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 500):
        return False, f"Molecular weight {mol_wt:.1f} Da not in expected range (150–500 Da) for sesquiterpenoids"
    
    # 4. Check aromaticity:
    aromatic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    if c_count > 0 and (aromatic_carbons / c_count) > 0.5:
        return False, "Too many aromatic carbons (more than 50%), which is not typical for a sesquiterpenoid structure"
    
    # 5. Extra checks for acyclic molecules.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rings == 0:
        # (a) Too flexible acyclic molecule – might be a fatty acid.
        if n_rotatable > 6:
            return False, f"Acyclic structure with high flexibility (rotatable bonds = {n_rotatable}), likely a fatty acid rather than a sesquiterpenoid"
        # (b) An almost linear carbon chain with a carboxylic acid group
        acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
        if mol.HasSubstructMatch(acid_smarts):
            longest_chain = _longest_carbon_chain(mol)
            if longest_chain >= (c_count - 2):
                return False, f"Acyclic fatty acid: longest carbon chain length ({longest_chain}) is almost the entire carbon count ({c_count}) in a molecule bearing a carboxylic acid"
    
    # 6. Count alkyl methyl groups.
    # We define an alkyl methyl group as a carbon atom with atomic number 6 that is:
    #  - bonded to exactly one other atom,
    #  - has exactly three hydrogens,
    #  - and that one neighbor is also a carbon (not O or S).
    alkyl_methyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6:
                alkyl_methyl_count += 1
    if alkyl_methyl_count < 3:
        return False, f"Too few alkyl methyl groups ({alkyl_methyl_count}) – sesquiterpenoids typically have at least 3 CH3 groups attached to carbons"
    
    return True, f"Molecule has {c_count} carbons, MW {mol_wt:.1f} Da, {alkyl_methyl_count} alkyl methyl groups, and passes allowed element and structural checks, consistent with a sesquiterpene-derived structure"

# Example usage (for testing purposes)
if __name__ == '__main__':
    tests = {
        "indicanone": "CC(=C)[C@@H]1CCC(CO)=C2CC(=O)C(C)=C2C1",
        "juvenile hormone I acid": "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",
        "(-)-GR24": "O1C(/C(/[C@H]2[C@@H]1C3=C(C2)C=CC=C3)=C/O[C@H]4OC(C(=C4)C)=O)=O",
        "3,7,11-trimethyl-dodecanoic acid": "OC(=O)CC(CCCC(CCCC(C)C)C)",
        "1,10,11,12-Guaianetetrol": "S(=O)(=O)(OCC(O)([C@H]1C[C@@H]2C(O)(CC[C@@H]2C)[C@@](CC1)(O)C)C)O",
        "geranyl acetone (false positive expected)": "CC(=O)CC\\C=C(/C)CCC=C(C)C",
        "Kadsurenone (false positive expected)": "COc1ccc(cc1OC)[C@H]1OC2=CC(=O)C(CC=C)=C[C@]2(OC)[C@@H]1C",
    }
    for name, smi in tests.items():
        result, reason = is_sesquiterpenoid(smi)
        print(f"{name}: {result} -- {reason}")