"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: CHEBI:49613 peptide antibiotic
A chemically diverse class of peptides that exhibit antimicrobial properties.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pickle

# Load pre-trained model
model = pickle.load(open("peptide_antibiotic_model.pkl", "rb"))

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string and a pre-trained machine learning model.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide backbone
    peptide_pattern = Chem.MolFromSmarts("[N;X3]([C;X4]([N;X3])=O)[C;X4][C;X3]")
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide backbone found"
    
    # Calculate molecular descriptors
    calc = MoleculeDescriptors.MoleculeDescriptors([mol])
    descriptors = calc.CalcDescriptors()

    # Predict using pre-trained model
    prediction = model.predict([descriptors])[0]
    
    if prediction == 1:
        return True, "Structural features and properties consistent with peptide antibiotic"
    else:
        return False, "Structural features and properties not consistent with peptide antibiotic"

__metadata__ = {
    'chemical_class': {'id': 'CHEBI:49613', 
                       'name': 'peptide antibiotic', 
                       'definition': 'A chemically diverse class of peptides that exhibit antimicrobial properties.', 
                       'parents': ['CHEBI:35827', 'CHEBI:36336']},
    'config': { 'llm_model_name': 'lbl/claude-sonnet',
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 2,
    'num_true_negatives': 182399,
    'num_false_negatives': 33,
    'num_negatives': None, 
    'precision': 0.9867549668874173,
    'recall': 0.8197674418604652,
    'f1': 0.9004329004329005,
    'accuracy': 0.9998522576737073}