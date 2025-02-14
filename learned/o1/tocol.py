"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: CHEBI:26416 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton that is substituted at position 2
    by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define chromanol core SMARTS pattern with atom mapping
    chromanol_smarts = """
    [O:1]1CC[C@:2]2(c3cc(O)ccc3OC2C1)  # Chroman-6-ol core with mapped atom at position 2
    """
    chromanol_pattern = Chem.MolFromSmarts(chromanol_smarts)
    if chromanol_pattern is None:
        return False, "Invalid chromanol SMARTS pattern"

    # Try to find chromanol core in the molecule
    chromanol_matches = mol.GetSubstructMatches(chromanol_pattern)
    if not chromanol_matches:
        return False, "No chromanol core (chroman-6-ol skeleton) found"

    # For each match, check substitution at position 2
    for match in chromanol_matches:
        # Map atom indices from SMARTS to molecule
        atom_map = {}
        for smarts_idx, mol_idx in enumerate(match):
            atom = chromanol_pattern.GetAtomWithIdx(smarts_idx)
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                atom_map[map_num] = mol_idx

        # Get position 2 atom (mapped as :2 in SMARTS)
        position_2_idx = atom_map.get(2)
        if position_2_idx is None:
            continue  # Try next match

        # Get substituents at position 2
        position_2_atom = mol.GetAtomWithIdx(position_2_idx)
        substituents = [nbr for nbr in position_2_atom.GetNeighbors() if nbr.GetIdx() not in match]
        if not substituents:
            continue  # No substitution at position 2

        # Assume the substituent is the side chain
        side_chain_atom = substituents[0]

        # Use BFS to get all atoms in the side chain
        visited = set()
        to_visit = [side_chain_atom.GetIdx()]
        side_chain_atoms = set()
        while to_visit:
            current_idx = to_visit.pop()
            if current_idx not in visited:
                visited.add(current_idx)
                side_chain_atoms.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)
                for nbr in current_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in visited and nbr_idx not in match:
                        to_visit.append(nbr_idx)

        # Check that side chain consists of only carbon and hydrogen atoms
        side_chain_elements = set(mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in side_chain_atoms)
        if not side_chain_elements.issubset({6, 1}):
            continue  # Side chain contains other atoms

        # Count number of carbons in side chain
        carbon_count = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count != 15:
            continue  # Side chain does not have 15 carbons (three isoprenoid units)

        # Identify isoprenoid pattern (three units of C5)
        isoprene_smarts = "[C;D1]([C;D2]=[C;D2])[C;D2]([C;D2])"
        isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
        side_chain_mol = Chem.PathToSubmol(mol, list(side_chain_atoms))
        isoprene_matches = side_chain_mol.GetSubstructMatches(isoprene_pattern)
        if len(isoprene_matches) not in [0, 3]:
            continue  # Side chain does not have proper isoprenoid units

        # Count number of double bonds in side chain
        double_bond_count = 0
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in side_chain_atoms and end_idx in side_chain_atoms:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_bond_count += 1

        if double_bond_count not in [0, 3]:
            continue  # Side chain is not saturated or triply-unsaturated

        return True, "Molecule is a tocol with chromanol core and appropriate side chain"

    # If no matches result in a tocol, return False
    return False, "Chromanol core or side chain does not meet criteria"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26416',
        'name': 'tocol',
        'definition': 'A chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}