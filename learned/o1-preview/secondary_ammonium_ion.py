"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:35277 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is an organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found_secondary_ammonium = False

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check if the nitrogen has a +1 formal charge
            if atom.GetFormalCharge() == 1:
                # Get neighbors of the nitrogen atom
                neighbors = atom.GetNeighbors()
                carbon_count = 0
                hydrogen_count = 0
                other_neighbors = 0
                valid_carbons = True

                for neighbor in neighbors:
                    atomic_num = neighbor.GetAtomicNum()
                    if atomic_num == 6:
                        carbon_count += 1
                        # Check if the carbon is part of a carbonyl group (exclude amides)
                        for bond in neighbor.GetBonds():
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                bonded_atom = bond.GetOtherAtom(neighbor)
                                if bonded_atom.GetAtomicNum() == 8:
                                    valid_carbons = False  # Carbon is carbonyl carbon
                    elif atomic_num == 1:
                        hydrogen_count += 1
                    else:
                        other_neighbors +=1  # Nitrogen bonded to other heteroatom

                # Check if nitrogen is bonded to exactly two carbons and two hydrogens
                if (carbon_count == 2 and hydrogen_count == 2 and other_neighbors == 0 and valid_carbons):
                    found_secondary_ammonium = True
                    return True, "Contains secondary ammonium ion group (protonated secondary amine)"
    
    if not found_secondary_ammonium:
        return False, "No secondary ammonium ion group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35277',
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['CHEBI:35274']
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
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}