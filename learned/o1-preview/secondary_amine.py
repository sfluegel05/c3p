"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary Amine
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is derived from ammonia where two hydrogen atoms are replaced by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            is_amine_nitrogen = True  # Assume it is an amine nitrogen

            # Exclude nitrogens in amides, nitro groups, nitriles, etc.
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    # Exclude amides (N-C=O) and nitro groups (N(=O)-O)
                    is_amine_nitrogen = False
                if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                    # Exclude nitriles (C#N)
                    is_amine_nitrogen = False
                if neighbor.GetAtomicNum() == 8 and neighbor.IsInRing():
                    # Exclude lactams (cyclic amides)
                    is_amine_nitrogen = False
                if neighbor.GetAtomicNum() == 16:
                    # Exclude sulfonamides (N-S(=O)2)
                    is_amine_nitrogen = False

            if not is_amine_nitrogen:
                continue  # Skip to the next nitrogen atom

            # Count the number of hydrogen atoms on nitrogen
            num_h = atom.GetTotalNumHs()

            # Get the formal charge
            formal_charge = atom.GetFormalCharge()

            # Get the degree (number of directly bonded atoms)
            degree = atom.GetDegree()

            # Get the list of neighbor atoms
            neighbors = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]

            # Check if nitrogen is connected to exactly two carbon atoms
            num_carbons = neighbors.count(6)
            if num_carbons != 2:
                continue  # Not a secondary amine

            # Check for secondary amine conditions
            # Neutral secondary amine: degree 3 (2 carbons + 1 hydrogen), formal charge 0
            if formal_charge == 0 and num_h == 1 and degree == 3:
                return True, "Molecule contains a secondary amine group"

            # Protonated secondary amine: degree 4 (2 carbons + 2 hydrogens), formal charge +1
            elif formal_charge == 1 and num_h == 2 and degree == 4:
                return True, "Molecule contains a protonated secondary amine group"

    return False, "No secondary amine group found"

__metadata__ = {'chemical_class': {'id': None,
                                   'name': 'secondary amine',
                                   'definition': 'A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.',
                                   'parents': []},
                'config': {'llm_model_name': 'lbl/claude-sonnet',
                           'f1_threshold': 0.8,
                           'max_attempts': 5,
                           'max_positive_instances': None,
                           'max_positive_to_test': None,
                           'max_negative_to_test': None,
                           'max_positive_in_prompt': 50,
                           'max_negative_in_prompt': 20,
                           'max_instances_in_prompt': 100,
                           'test_proportion': 0.1},
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
                'accuracy': None}