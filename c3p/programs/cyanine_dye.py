"""
Classifies: CHEBI:37960 cyanine dye
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cyanine_dye(smiles: str):
    """
    Determines if a molecule is a cyanine dye.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanine dye, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of conjugated chain R2N[CH=CH]nCH=N(+)R2 <-> R2N(+)=CH[CH=CH]nNR2
    pattern = Chem.MolFromSmarts("N([#6])[#6]=[#6][#6]=[#6][#6]=[#6][#6]=[#6]N(+)")  # General pattern for cyanine dyes
    if mol.HasSubstructMatch(pattern):
        return True, "Contains the conjugated chain characteristic of cyanine dyes"

    return False, "Does not contain the conjugated chain characteristic of cyanine dyes"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37960',
                          'name': 'cyanine dye',
                          'definition': 'Cyanine dyes are synthetic dyes with '
                                        'the general formula '
                                        'R2N[CH=CH]nCH=N(+)R2    <-> '
                                        'R2N(+)=CH[CH=CH]nNR2 (n is a small '
                                        'number) in which the nitrogen and '
                                        'part of the conjugated chain usually '
                                        'form part of a heterocyclic system, '
                                        'such as imidazole, pyridine, pyrrole, '
                                        'quinoline and thiazole.',
                          'parents': ['CHEBI:35352']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': 'Python argument types in\n'
             '    Mol.HasSubstructMatch(Mol, NoneType)\n'
             'did not match C++ signature:\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'RDKit::SubstructMatchParameters params=True)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'RDKit::SubstructMatchParameters params)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}