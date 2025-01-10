"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are esters formed from long-chain fatty acids and long-chain alcohols,
    typically with unbranched carbon chains of 14 carbons or more on both sides of the ester linkage.
    They do not have additional functional groups or rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester functional groups
    ester_smarts = "[#6][C](=O)[O][#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_esters = len(ester_matches)
    if num_esters != 1:
        return False, f"Found {num_esters} ester groups, expected exactly 1"

    # Check for other functional groups
    # Define SMARTS patterns for common functional groups
    func_groups = {
        'acid': '[CX3](=O)[OX1H0-,OX2H1]',
        'aldehyde': '[CX3H1](=O)[#6]',
        'ketone': '[#6][CX3](=O)[#6]',
        'alcohol': '[#6][OX2H]',
        'amine': '[NX3;H2,H1;!$(NC=O)]',
        'amide': '[NX3][CX3](=O)[#6]',
        'ether': '[CX4][OX2][CX4]',
        'nitrile': '[CX2]#N',
        'halogen': '[F,Cl,Br,I]',
        'sulfide': '[#16]',
        'phosphate': '[P]',
        'sulfonamide': 'S(=O)(=O)N',
        'nitro': '[$([NX3](=O)=O)]'
    }
    for name, pattern in func_groups.items():
        smarts = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(smarts)
        if matches:
            return False, f"Contains {name} functional group"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structures"

    # Get the ester atoms
    ester_match = ester_matches[0]
    carbonyl_c_idx = ester_match[1]
    single_o_idx = ester_match[2]
    alkyl_c_idx = ester_match[3]

    # Function to count carbons in a chain recursively
    def count_chain_carbons(atom_idx, visited=None):
        if visited is None:
            visited = set()
        if atom_idx in visited:
            return 0
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return 0  # Only count carbons
        count = 1
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            # Skip backtracking to ester oxygen or carbon
            if neighbor_idx in visited:
                continue
            neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
            if neighbor_atom.GetAtomicNum() in (6, 1):  # Carbon or hydrogen
                count += count_chain_carbons(neighbor_idx, visited)
            else:
                return 0  # Branching to heteroatom, not acceptable
        return count

    # Count carbons on acyl side (excluding ester carbonyl carbon)
    acyl_chain_carbons = 0
    for neighbor in mol.GetAtomWithIdx(carbonyl_c_idx).GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx != ester_match[0]:
            acyl_chain_carbons = count_chain_carbons(neighbor_idx)
            break

    # Count carbons on alkoxy side
    alkoxy_chain_carbons = count_chain_carbons(alkyl_c_idx)

    # Check for branching in acyl chain
    acyl_chain_atoms = Chem.PathToSubmol(mol, mol.GetSubstructMatch(Chem.MolFromSmarts(f"[C:1][C:2](=O)[O][C]"), uniquify=False))
    if acyl_chain_atoms and acyl_chain_atoms.GetNumAtoms() != acyl_chain_carbons + 2:
        return False, "Acyl chain is branched"

    # Check for branching in alkoxy chain
    alkoxy_chain_atoms = Chem.PathToSubmol(mol, mol.GetSubstructMatch(Chem.MolFromSmarts(f"[O][C:1][C]"), uniquify=False))
    if alkoxy_chain_atoms and alkoxy_chain_atoms.GetNumAtoms() != alkoxy_chain_carbons + 1:
        return False, "Alkoxy chain is branched"

    if acyl_chain_carbons < 14:
        return False, f"Acyl chain too short ({acyl_chain_carbons} carbons)"
    if alkoxy_chain_carbons < 14:
        return False, f"Alkoxy chain too short ({alkoxy_chain_carbons} carbons)"

    return True, f"Valid wax: ester with unbranched acyl chain ({acyl_chain_carbons} carbons) and unbranched alkoxy chain ({alkoxy_chain_carbons} carbons)"

__metadata__ = {
    'chemical_class': {
        'name': 'wax',
        'definition': 'A chemical substance that is an organic compound or mixture of compounds that is composed of long-chain molecules and is malleable at ambient temperatures.'
    }
}