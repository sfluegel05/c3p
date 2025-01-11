"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26195 phytosterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol, occurring in plants, with variations in carbon side chains and/or presence/absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid nucleus pattern (tetracyclic ring system with a hydroxyl group at the 3-position)
    steroid_nucleus_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C")
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "No steroid nucleus found (tetracyclic ring system with a hydroxyl group at the 3-position)"

    # Check for a side chain at the 17-position (typical of sterols)
    side_chain_pattern = Chem.MolFromSmarts("[C@@]12CC[C@@H]3[C@@]4(CC[C@H](C4)O)C[C@@H]1CC[C@@]2(C)CC")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No side chain found at the 17-position"

    # Check for variations in the side chain (e.g., double bonds, additional carbons)
    # Phytosterols often have longer or unsaturated side chains compared to cholesterol
    side_chain_atoms = mol.GetSubstructMatch(side_chain_pattern)
    side_chain_mol = Chem.PathToSubmol(mol, side_chain_atoms)
    if side_chain_mol is None:
        return False, "Unable to extract side chain for analysis"

    # Count the number of carbons in the side chain
    side_chain_c_count = sum(1 for atom in side_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if side_chain_c_count < 8:
        return False, "Side chain too short for a typical phytosterol"

    # Check for double bonds in the side chain (common in phytosterols)
    side_chain_double_bonds = sum(1 for bond in side_chain_mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if side_chain_double_bonds == 0:
        return False, "No double bonds found in the side chain (common in phytosterols)"

    # Check molecular weight (phytosterols typically have a higher molecular weight than cholesterol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a typical phytosterol"

    return True, "Contains a steroid nucleus with a hydroxyl group at the 3-position and a varied side chain at the 17-position"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26195',
        'name': 'phytosterol',
        'definition': 'Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.',
        'parents': ['CHEBI:26195', 'CHEBI:26195']
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}