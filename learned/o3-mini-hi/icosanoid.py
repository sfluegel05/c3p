"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation
of the three C20 essential fatty acids (EFAs) – icosapentaenoic acid (EPA), arachidonic acid (AA)
and dihomo-gamma-linolenic acid (DGLA).

Improved heuristic criteria:
  • SMILES must be valid.
  • Molecular weight between 250 and 900 Da.
  • At least 3 oxygen atoms.
  • At least one and at most 4 non‐aromatic C=C bonds.
  • The longest continuous carbon chain (via DFS on carbon atoms) between 15 and 22 atoms.
  • The overall ring count must be ≤ 2.
  • If the molecule does not contain phosphorus then the presence of any nitrogen atoms suggests a non‐icosanoid.
  
Note: This heuristic is an approximation and may mis‐classify borderline cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the longest continuous chain (i.e. simple path) of carbon atoms.
    We restrict the graph to carbon (atomic number 6) and perform DFS from each node.
    """
    # Build a graph for carbon atoms (by index)
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx = atom.GetIdx()
            carbon_graph[idx] = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    carbon_graph[idx].append(nbr.GetIdx())
    
    max_length = 0
    def dfs(node, visited):
        current_max = 1
        for neighbor in carbon_graph.get(node, []):
            if neighbor not in visited:
                new_length = 1 + dfs(neighbor, visited | {neighbor})
                current_max = max(current_max, new_length)
        return current_max

    for node in carbon_graph:
        chain_length = dfs(node, {node})
        if chain_length > max_length:
            max_length = chain_length
    return max_length

def is_icosanoid(smiles: str):
    """
    Determines whether a given molecule qualifies as an icosanoid based on its SMILES string.
    
    Heuristic filters applied:
      - Valid SMILES.
      - Molecular weight between 250 and 900 Da.
      - Contains at least 3 oxygen atoms.
      - Contains at least 1 and no more than 4 non‐aromatic C=C double bonds.
      - The longest continuous carbon chain is between 15 and 22 atoms.
      - The total ring count is ≤ 2.
      - For molecules that do not contain phosphorus, the presence of any nitrogen atoms 
        indicates the molecule is not an icosanoid.
    
    Additionally, we check whether a free carboxylic acid is present (many icosanoids are acids, 
    though conjugated/esterified examples exist).
    
    Returns:
        (bool, str): True if the molecule meets the heuristic criteria with an explanation,
                     False otherwise with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count total carbons and oxygens.
    total_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the expected range (250-900 Da)."
    
    # Oxygen count should be at least 3.
    if oxygen_count < 3:
        return False, f"Insufficient oxygen count ({oxygen_count}); typical icosanoids have ≥3 oxygen atoms."
    
    # Count non-aromatic C=C double bonds.
    alkene_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and (not bond.GetIsAromatic()):
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                alkene_count += 1
    if alkene_count < 1:
        return False, "Too few non-aromatic C=C double bonds; icosanoids are typically unsaturated."
    if alkene_count > 4:
        return False, f"Too many non-aromatic C=C double bonds ({alkene_count}); not typical for icosanoids."
    
    # Compute the longest continuous carbon chain.
    longest_chain = longest_carbon_chain(mol)
    if longest_chain < 15 or longest_chain > 22:
        return False, (f"Longest continuous carbon chain ({longest_chain}) is not in the expected range (15-22).")
    
    # Count the total ring count.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings > 2:
        return False, f"Ring count ({n_rings}) is too high; icosanoids typically have ≤2 rings."
    
    # Check heteroatom composition.
    n_N = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    n_P = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    # If no phosphorus is present and there is any nitrogen, then we suspect a different class.
    if n_P == 0 and n_N > 0:
        return False, f"Non-conjugated molecule contains nitrogen (N count: {n_N}); not typical for icosanoids."
    
    # Check for free carboxylic acid.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if mol.HasSubstructMatch(acid_pattern):
        acid_note = ""
    else:
        acid_note = " (No free carboxylic acid detected; molecule may be fully esterified/conjugated.)"
    
    reason = (f"Meets heuristic criteria for an icosanoid: (Total C={total_carbon}, O={oxygen_count}, "
              f"MW={mol_wt:.1f} Da, non-aromatic C=C count={alkene_count}, Longest C-chain={longest_chain}, "
              f"Rings={n_rings}){acid_note}")
    
    return True, reason

# Example usage and test cases:
if __name__ == "__main__":
    test_smiles = [
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4 (true positive)
        "CCCC[C@H](O)\\C=C\\[C@H]([C@@H](C\\C=C/CCCC(O)=O)C=O)C(C)=O",  # levuglandin D2 (true positive)
        "O[C@H]1[C@@H]([C@H]([C@@H](O)C1)C/C=C\\CCCC(OC(C)C)=O)CC[C@@H](O)CCOC2=CC=CC=C2",  # 17-phenoxy trinor PGF2alpha ester (true positive)
        "CCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O",  # (8R)-hydroxy-(11S,12S)-epoxyicosa… acid (true positive)
        "[H]C(CC)=C([H])CCCCCCCCCCCCCCCC(O)=O",  # 17-icosenoic acid (false negative: oxygen count too low)
        "CCCC[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",  # 15(S)-HPETE (true positive)
        # A few false positive examples based on heteroatom/ring criteria:
        "O=C(N[C@@H](C(O)(C)C)C)[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C",  # Contains nitrogen not allowed
        "COc1ccc(cc1OC)[C@H]1OC2=CC(=O)C(CC=C)=C[C@]2(OC)[C@@H]1C",  # Kadsurenone, rings likely >2
    ]
    
    for smi in test_smiles:
        result, explanation = is_icosanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")