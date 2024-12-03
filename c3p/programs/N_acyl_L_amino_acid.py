"""
Classifies: CHEBI:21644 N-acyl-L-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an amino acid backbone (N-C-C(=O)-O)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[C@@H](C(=O)O)N')):
        return False, "No L-amino acid backbone found"

    # Check for the presence of an acyl group attached to the nitrogen
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('NC(=O)')):
        return False, "No acyl group attached to the nitrogen"

    # Check for L-configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not any(center[1] == 'S' for center in chiral_centers):
        return False, "No L-configuration found"

    return True, "Molecule is an N-acyl-L-amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:21644',
                          'name': 'N-acyl-L-amino acid',
                          'definition': 'Any N-acylamino acid having '
                                        'L-configuration.',
                          'parents': ['CHEBI:51569']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 5,
    'precision': 0.8823529411764706,
    'recall': 0.75,
    'f1': 0.8108108108108107,
    'accuracy': None}