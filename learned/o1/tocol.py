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

    # Define chromanol core SMARTS pattern
    chromanol_smarts = "c1cc(O)ccc1C2CCCO2"  # Chromanol core (chroman-6-ol skeleton)
    chromanol_pattern = Chem.MolFromSmarts(chromanol_smarts)
    if chromanol_pattern is None:
        return False, "Invalid chromanol SMARTS pattern"

    # Find chromanol core matches
    matches = mol.GetSubstructMatches(chromanol_pattern)
    if not matches:
        return False, "No chromanol core (chroman-6-ol skeleton) found"

    # For each match, check substitution at position 2
    for match in matches:
        # Get the atoms in the match
        match_atoms = list(match)
        # Identify the oxygen atom in the chromanol core
        oxygen_atom_idx = None
        for idx in match_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8 and atom.IsInRing():
                oxygen_atom_idx = idx
                break
        if oxygen_atom_idx is None:
            continue  # No oxygen atom found in ring
        # Get the carbon atom at position 2 (adjacent to oxygen in ring)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_atom_idx)
        ring_carbon_idxs = [nbr.GetIdx() for nbr in oxygen_atom.GetNeighbors() if nbr.GetIdx() in match_atoms and nbr.GetAtomicNum()==6]
        if len(ring_carbon_idxs) < 2:
            continue  # Not enough ring carbons found
        # Check for substituents on ring carbons
        position_2_atom = None
        for carbon_idx in ring_carbon_idxs:
            carbon_atom = mol.GetAtomWithIdx(carbon_idx)
            substituents = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetIdx() not in match_atoms]
            if substituents:
                position_2_atom = carbon_atom
                break
        if position_2_atom is None:
            continue  # No substitution at position 2
        # Get the side chain starting from the substituent
        side_chain_atoms = set()
        to_visit = [nbr.GetIdx() for nbr in position_2_atom.GetNeighbors() if nbr.GetIdx() not in match_atoms]
        while to_visit:
            current_idx = to_visit.pop()
            if current_idx not in side_chain_atoms:
                side_chain_atoms.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)
                for nbr in current_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in side_chain_atoms and nbr_idx not in match_atoms:
                        to_visit.append(nbr_idx)
        # Check that side chain consists of only carbon and hydrogen atoms
        side_chain_elements = set(mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in side_chain_atoms)
        if not side_chain_elements.issubset({6, 1}):
            continue  # Side chain contains other atoms
        # Count number of carbons in side chain
        carbon_count = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count != 15:
            continue  # Side chain does not have 15 carbons (three isoprenoid units)
        # Count number of double bonds in side chain
        double_bond_count = 0
        for idx in side_chain_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                neighbor_idx = bond.GetOtherAtomIdx(idx)
                if neighbor_idx in side_chain_atoms and idx < neighbor_idx:
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bond_count += 1
        if double_bond_count not in [0, 3]:
            continue  # Side chain is not saturated or triply unsaturated
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}