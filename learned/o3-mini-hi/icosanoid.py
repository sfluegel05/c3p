"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation
of the three C20 essential fatty acids (EFAs) – icosapentaenoic acid (EPA), arachidonic acid (AA)
and dihomo-gamma-linolenic acid (DGLA).

Heuristic criteria applied (improved):
  • The SMILES must be valid.
  • Molecular weight between 250 and 900 Da.
  • At least 3 oxygen atoms (to help avoid simple fatty acids).
  • Must have at least one non‐aromatic C=C double bond.
  • The longest continuous carbon chain (via DFS on only carbon atoms) must be between 15 and 22 carbons.
  • For molecules that are fully linear (i.e. no rings and longest chain equals total carbon count),
    at least 2 non‐aromatic C=C bonds are required.
    
Note: This heuristic – though improved relative to the previous version – is still an approximation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Compute the longest continuous chain (simple path) consisting solely of carbon atoms.
    We build a graph for all carbon atoms (atomic number 6) and perform DFS from each node.
    """
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
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
        max_length = max(max_length, chain_length)
    return max_length

def is_icosanoid(smiles: str):
    """
    Determines whether a molecule qualifies as an icosanoid based on its SMILES string.
    
    Improved heuristic filters:
      - Valid SMILES.
      - Molecular weight between 250 and 900 Da.
      - Must contain at least 3 oxygen atoms (helps dismiss simple fatty acids).
      - Contains at least one non-aromatic C=C double bond.
      - The longest continuous carbon chain must be between 15 and 22 carbons.
      - For fully linear molecules (no rings and longest chain equals total carbon count),
        require ≥2 non-aromatic C=C bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the icosanoid heuristic criteria; False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Count total carbons and oxygens.
    total_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the expected range (250-900 Da)."
    
    # Check for sufficient oxygenation; simple fatty acids usually have too few oxygens.
    if oxygen_count < 3:
        return False, f"Insufficient oxygen count ({oxygen_count}); typical icosanoids have ≥3 oxygen atoms."
    
    # Count non-aromatic C=C double bonds.
    alkene_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                alkene_count += 1
    if alkene_count < 1:
        return False, "Too few non-aromatic C=C double bonds; icosanoids are typically unsaturated."
    
    # Compute the longest continuous carbon chain.
    longest_chain = longest_carbon_chain(mol)
    if longest_chain < 15 or longest_chain > 22:
        return False, (f"Longest continuous carbon chain ({longest_chain}) is not in the expected range "
                       "(15-22).")
    
    # Get ring count.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    
    # If molecule is fully linear then require at least 2 non-aromatic C=C bonds to avoid simple fatty acids.
    if n_rings == 0 and longest_chain == total_carbon:
        if alkene_count < 2:
            return False, ("Molecule appears fully linear (no rings, chain equals total carbons) with only one C=C; "
                           "likely a simple fatty acid rather than an oxidized icosanoid.")
    
    # Note on acid group: Many icosanoids are free acids, yet some conjugates/esters are also known.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        acid_note = " (No free carboxylic acid detected; molecule may be fully esterified/conjugated.)"
    else:
        acid_note = ""
    
    reason = (f"Meets heuristic criteria for an icosanoid: (Total C={total_carbon}, O={oxygen_count}, "
              f"MW={mol_wt:.1f} Da, non-aromatic C=C count={alkene_count}, Longest C-chain={longest_chain}, "
              f"Rings={n_rings}){acid_note}")
    return True, reason

# Example usage:
if __name__ == "__main__":
    # Test a set of example SMILES from the provided list.
    smiles_examples = [
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCC[C@H](O)\\C=C\\[C@H]([C@@H](C\\C=C/CCCC(O)=O)C=O)C(C)=O",  # levuglandin D2
        "O[C@H]1[C@@H]([C@H]([C@@H](O)C1)C/C=C\\CCCC(OC(C)C)=O)CC[C@@H](O)CCOC2=CC=CC=C2",  # 17-phenoxy trinor Prostaglandin F2alpha isopropyl ester
        "CCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O",  # (8R)-hydroxy-(11S,12S)-epoxyicosa-(5Z,9E,14Z)-trienoic acid
        "[H]C(CC)=C([H])CCCCCCCCCCCCCCCC(O)=O",  # 17-icosenoic acid (expected to be rejected)
        "CCCC[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",  # 15(S)-HPETE
    ]
    
    for smi in smiles_examples:
        result, explanation = is_icosanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")