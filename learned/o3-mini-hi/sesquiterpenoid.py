"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: Sesquiterpenoid
Definition:
  Any terpenoid derived from a sesquiterpene (C15 skeleton) even if the skeleton
  has been rearranged or modified (e.g., by removal of one or more skeletal atoms).
  
Heuristic improvements (this version):
  - Allowed elements: only hydrogen, carbon, oxygen, and sulfur.
  - Carbon count is required to be between 14 and 21.
  - Molecular weight must be between 150 and 500 Da.
  - The degree of unsaturation (DBE) is computed and is required to be between 2 and 7.
  - The fraction of aromatic carbons must be ≤50%.
  - For acyclic molecules: if the structure is very flexible (>6 rotatable bonds)
    or has an almost linear carbon chain with a carboxylic acid group, then it is likely a fatty acid.
  - Additionally, we require at least 2 “alkyl” methyl groups (CH3 groups bonded to carbon)
    in recognition that rearrangements may sometimes lower the “expected” count.
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
    # Only consider carbon atoms
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
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

def _calc_dbe(mol, c_count: int):
    """
    Calculate the degree of unsaturation (Double Bond Equivalents, DBE) using the 
    formula: DBE = (2*C + 2 - H)/2.
    Only carbons and hydrogens are considered since O and S do not change DBE.
    """
    h_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            h_count += 1
        # include explicit hydrogens attached to heavy atoms too
        h_count += atom.GetTotalNumHs()
    return (2 * c_count + 2 - h_count) / 2

def is_sesquiterpenoid(smiles: str):
    """
    Determines whether a molecule is a sesquiterpenoid (or derived from a sesquiterpene)
    using heuristic rules.

    The following checks are used:
      1. Allowed elements: only H, C, O, and S.
      2. Carbon count must be roughly 14–21.
      3. Molecular weight must lie between 150 and 500 Da.
      4. Degree of unsaturation (DBE) must be between 2 and 7.
      5. The fraction of aromatic carbons must be ≤50%.
      6. If the molecule is acyclic and either very flexible (rotatable bonds >6) 
         or has an almost unbranched carbon chain with a carboxylic acid group,
         then it is considered likely a fatty acid rather than a sesquiterpenoid.
      7. The molecule should have at least 2 “alkyl” methyl groups (CH3 groups bonded to a carbon)
         since terpene-derived structures typically feature several methyl substituents.
         
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple where the boolean indicates if the molecule matches the heuristic and
                   the string gives a reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Allowed elements: only H, C, O, and S.
    allowed_atomic_nums = {1, 6, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()}, not allowed for a typical sesquiterpenoid"
    
    # 2. Carbon count check.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (14 <= c_count <= 21):
        return False, f"Molecule contains {c_count} carbons, not consistent with a typical sesquiterpene-derived structure (expected roughly 14–21 carbons)"
    
    # 3. Molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 500):
        return False, f"Molecular weight {mol_wt:.1f} Da not in expected range (150–500 Da) for sesquiterpenoids"
    
    # 4. Degree of unsaturation (DBE) check.
    dbe = _calc_dbe(mol, c_count)
    if not (2 <= dbe <= 7):
        return False, f"Degree of unsaturation (DBE = {dbe:.1f}) not in expected range (2–7) for sesquiterpenoids"
    
    # 5. Aromatic carbons fraction.
    aromatic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    if c_count > 0 and (aromatic_carbons / c_count) > 0.5:
        return False, "Too many aromatic carbons (more than 50%), which is not typical for a sesquiterpenoid structure"
    
    # 6. Acyclic checks: if no rings, then check for flexibility or linear acid chain.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rings == 0:
        if n_rotatable > 6:
            return False, f"Acyclic structure with high flexibility (rotatable bonds = {n_rotatable}), likely a fatty acid rather than a sesquiterpenoid"
        acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
        if mol.HasSubstructMatch(acid_smarts):
            longest_chain = _longest_carbon_chain(mol)
            if longest_chain >= (c_count - 2):
                return False, f"Acyclic fatty acid: longest carbon chain length ({longest_chain}) is almost the entire carbon count ({c_count}) in a molecule bearing a carboxylic acid"
    
    # 7. Count alkyl methyl groups.
    # Define an alkyl methyl group as a carbon with:
    #    - atomic number 6, degree 1, exactly three hydrogens,
    #    - and its sole neighbor is also carbon.
    alkyl_methyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6:
                alkyl_methyl_count += 1
    if alkyl_methyl_count < 2:
        return False, f"Too few alkyl methyl groups ({alkyl_methyl_count}) – sesquiterpenoids typically have at least 2 CH3 groups attached to carbons"
    
    return True, (f"Molecule has {c_count} carbons, MW {mol_wt:.1f} Da, DBE {dbe:.1f}, "
                  f"{alkyl_methyl_count} alkyl methyl groups, and passes allowed element and structural checks, "
                  "consistent with a sesquiterpene-derived structure")

# Example usage (for testing purposes)
if __name__ == '__main__':
    tests = {
        "indicanone": "CC(=C)[C@@H]1CCC(CO)=C2CC(=O)C(C)=C2C1",
        "juvenile hormone I acid": "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",
        "(-)-GR24": "O1C(/C(/[C@H]2[C@@H]1C3=C(C2)C=CC=C3)=C/O[C@H]4OC(C(=C4)C)=O)=O",
        "3,7,11-trimethyl-dodecanoic acid": "OC(=O)CC(CCCC(CCCC(C)C)C)",
        "1,10,11,12-Guaianetetrol": "S(=O)(=O)(OCC(O)([C@H]1C[C@@H]2C(O)(CC[C@@H]2C)[C@@](CC1)(O)C)C)O",
        "Trepibutone (false positive expected?)": "C(CCC(O)=O)(=O)C1=C(OCC)C=C(C(=C1)OCC)OCC",
    }
    for name, smi in tests.items():
        result, reason = is_sesquiterpenoid(smi)
        print(f"{name}: {result} -- {reason}")