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
  • (Total carbon count is noted, and if the molecule is completely linear with no rings, 
    it is assumed to be a simple fatty acid.)
  • Contains at least one oxygenated carbonyl group. Here the SMARTS "[CX3]=O" is used,
    unless the molecule contains phosphorus (which indicates possible conjugation).
  • Contains at least one non‐aromatic C=C (alkene) double bond.
  • The longest continuous carbon chain (computed via DFS) must be at least 15 carbons long.
    
Note: This heuristic – though improved relative to the previous version – is still an approximation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Compute the longest simple path (chain) consisting solely of carbon atoms.
    We first build a graph for carbon–atoms (atomic number 6) and then perform a DFS.
    """
    # Build a graph: for every carbon, list neighboring carbons.
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx = atom.GetIdx()
            carbon_graph[idx] = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    carbon_graph[idx].append(nbr.GetIdx())
                    
    # DFS function to find the longest simple path.
    max_length = 0
    def dfs(node, visited):
        current_max = 1
        for neighbor in carbon_graph.get(node, []):
            if neighbor not in visited:
                new_length = 1 + dfs(neighbor, visited | {neighbor})
                if new_length > current_max:
                    current_max = new_length
        return current_max

    for node in carbon_graph:
        chain_length = dfs(node, {node})
        if chain_length > max_length:
            max_length = chain_length
    return max_length

def is_icosanoid(smiles: str):
    """
    Determines whether a molecule qualifies as an icosanoid based on its SMILES string.
    
    Heuristic filters:
      - The SMILES must be valid.
      - Molecular weight between 250 and 900 Da.
      - Count total carbon atoms (with note: if the molecule is fully linear – no rings
        and the longest chain equals the total carbon count – it is likely a simple fatty acid).
      - Must contain at least one oxygenated carbonyl fragment. We use a general SMARTS "[CX3]=O"
        unless there is a phosphorus atom present (which may indicate conjugation).
      - Must have at least one non-aromatic C=C double bond.
      - The longest continuous carbon chain must be at least 15 carbons long.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the heuristic for icosanoids, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Total number of carbon atoms.
    total_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside expected 250-900 Da range."
    
    # Check for at least one oxygenated carbonyl fragment.
    # Using a more general SMARTS that matches any carbonyl group.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    # If the molecule contains a phosphorus atom, we assume the carbonyl might be masked (e.g. conjugated ester/phosphate)
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
    
    # Compute longest continuous carbon chain.
    longest_chain = longest_carbon_chain(mol)
    if longest_chain < 15:
        return False, f"Longest continuous carbon chain length ({longest_chain}) is too short to represent a fatty acid backbone."
    
    # Check for "simple fatty acid" pattern:
    # If there are no rings and the longest continuous chain equals the total carbon count,
    # then the molecule is likely just a linear fatty acid.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings == 0 and longest_chain == total_carbon:
        return False, "Molecule appears to be a simple fatty acid (fully linear with no rings)."
    
    # Optionally, look for a free carboxylic acid pattern.
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
    # A selection of example SMILES strings (from the provided list)
    smiles_examples = [
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCC[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",      # 15(S)-HPETE
        "O[C@H]1[C@@H]([C@H]([C@@H](O)C1)C/C=C\\CCCC(OC(C)C)=O)CC[C@@H](O)CCOC2=CC=CC=C2",  # 17-phenoxy trinor PGF2alpha isopropyl ester
        "[H]C(CC)=C([H])CCCCCCCCCCCCCCCC(O)=O",               # 17-icosenoic acid (icosanoid, though esterified/conjugated note may appear)
        "CCCCCC\\C=C\\CCCCCCCCC([O-])=O"                      # (10E)-octadecenoate, a simple linear fatty acid (should be rejected)
    ]
    
    for smi in smiles_examples:
        result, reason = is_icosanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")