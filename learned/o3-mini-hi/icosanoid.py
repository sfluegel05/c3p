"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation
of the three C20 essential fatty acids (EFAs) – icosapentaenoic acid (EPA), arachidonic acid (AA)
and dihomo-gamma-linolenic acid (DGLA).

Heuristic criteria applied:
  • The SMILES must be valid.
  • Molecular weight between 250 and 900 Da.
  • Requires at least one oxygenated carbonyl group ([CX3]=O) unless phosphorus is present.
  • Must have at least one non‐aromatic C=C double bond.
  • The longest continuous carbon chain (via DFS on only carbon atoms) must be at least 15 atoms.
  • For molecules that are fully linear (no rings and the longest chain equals total carbon count),
    we require at least 2 non‐aromatic C=C bonds (to distinguish polyunsaturated, oxidized species
    from simple fatty acids).
    
Note: This heuristic – though improved relative to the previous version – is still an approximation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Compute the longest continuous chain (simple path) consisting solely of carbon atoms.
    This is done by building a graph for all carbon atoms and doing a DFS on each.
    """
    carbon_graph = {}
    # Build a graph that maps indices of carbon atoms to their carbon neighbors.
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
    
    Heuristic filters:
      - Valid SMILES.
      - Molecular weight between 250 and 900 Da.
      - Contains at least one oxygenated carbonyl fragment ([CX3]=O) (unless phosphorus is present).
      - Has at least one non-aromatic C=C double bond.
      - The longest continuous carbon chain must be at least 15 carbons long.
      - For molecules that are fully linear (no rings and longest chain equals total carbon count),
        at least 2 non-aromatic C=C bonds are required to help distinguish oxidized fatty acids
        (icosanoids) from simple fatty acids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the heuristic for icosanoids, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count total number of carbon atoms.
    total_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the expected 250-900 Da range."
    
    # Check for at least one oxygenated carbonyl fragment ([CX3]=O).
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    contains_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if not contains_phosphorus and not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No oxygenated carbonyl fragment ([CX3]=O) found."
    
    # Count non-aromatic C=C bonds.
    alkene_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and
                not bond.GetIsAromatic()):
                alkene_count += 1
    if alkene_count < 1:
        return False, "Too few non-aromatic C=C double bonds; icosanoids are typically unsaturated."
    
    # Compute the longest continuous carbon chain.
    longest_chain = longest_carbon_chain(mol)
    if longest_chain < 15:
        return False, f"Longest continuous carbon chain length ({longest_chain}) is too short."
    
    # Get information about rings.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    
    # If molecule is fully linear (i.e. no rings and the longest chain equals total carbon count),
    # then require at least 2 non-aromatic C=C bonds to avoid classifying simple monounsaturated fatty acids.
    if n_rings == 0 and longest_chain == total_carbon:
        if alkene_count < 2:
            return False, ("Molecule appears to be fully linear with only one C=C; "
                           "likely a simple fatty acid rather than an oxidized icosanoid.")
    
    # Optionally, check and note if a free carboxylic acid group is present.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        acid_note = " (No free carboxylic acid detected; molecule may be fully esterified/conjugated.)"
    else:
        acid_note = ""
    
    reason = (f"Meets heuristic criteria for an icosanoid: "
              f"(Total C={total_carbon}, MW={mol_wt:.1f} Da, non-aromatic C=C count={alkene_count}, "
              f"Longest C-chain={longest_chain}, Rings={n_rings}){acid_note}")
    return True, reason

# Example usage:
if __name__ == "__main__":
    # A few sample SMILES from the provided list.
    smiles_examples = [
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCC[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",      # 15(S)-HPETE (was previously missed because fully linear; now requires >=2 alkenes)
        "CCCC\\C=C/C[C@@H](O)\\C=C\\C=C/C=C/[C@@H](O)CCCC(O)=O", # Delta(6)-trans,Delta(8)-cis-leukotriene B4 (fully linear, but with multiple C=C)
        "[H]C(CC)=C([H])CCCCCCCCCCCCCCCC(O)=O",               # 17-icosenoic acid (simple fatty acid – should be rejected)
        "CCCC[C@H](O)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",      # 11(R)-HETE (fully linear with multiple C=C)
    ]
    
    for smi in smiles_examples:
        result, reason = is_icosanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")