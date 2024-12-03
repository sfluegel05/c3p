"""
Classifies: CHEBI:39201 phenazines
"""
from rdkit import Chem

def is_phenazines(smiles: str):
    """
    Determines if a molecule is a phenazine or its derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenazine or derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least three rings
    if len(rings.AtomRings()) < 3:
        return False, "Less than three rings found"

    # Find the phenazine core structure (two benzene rings fused via two nitrogen atoms)
    phenazine_core = None
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetSymbol() == 'C' or atom.GetSymbol() == 'N' for atom in atoms):
                if sum(atom.GetSymbol() == 'N' for atom in atoms) == 2:
                    phenazine_core = ring
                    break

    if phenazine_core is None:
        return False, "No phenazine core structure found"

    # Check if the core structure is aromatic
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in phenazine_core):
        return False, "Phenazine core structure is not aromatic"

    return True, "Phenazine or its derivative found"

# Example usage
smiles = "O=C(N)C=1C2=NC=3C(=C(O)C=CC3)N=C2C=CC1"
print(is_phenazines(smiles))  # Expected output: (True, "Phenazine or its derivative found")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:39201',
                          'name': 'phenazines',
                          'definition': 'Any  organonitrogen heterocyclic '
                                        'compound based on a phenazine '
                                        'skeleton and derivatives.',
                          'parents': ['CHEBI:26979', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Phenazine or its derivative found')\n",
    'num_true_positives': 16,
    'num_false_positives': 2,
    'num_true_negatives': 17,
    'num_false_negatives': 3,
    'precision': 0.8888888888888888,
    'recall': 0.8421052631578947,
    'f1': 0.8648648648648649,
    'accuracy': None}