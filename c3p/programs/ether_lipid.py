"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for ether linkage (carbon-oxygen-carbon)
    ether_pattern = Chem.MolFromSmarts("COC")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No ether linkage found"

    # Check for alkyl chain attached to glycerol via ether linkage
    alkyl_chain_pattern = Chem.MolFromSmarts("COC[*]")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No alkyl chain attached to glycerol via ether linkage found"

    # Check for the presence of a phosphocholine or phosphoethanolamine group
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[NH3+]")
    if not (mol.HasSubstructMatch(phosphocholine_pattern) or mol.HasSubstructMatch(phosphoethanolamine_pattern)):
        return False, "No phosphocholine or phosphoethanolamine group found"

    return True, "Molecule is an ether lipid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64611',
                          'name': 'ether lipid',
                          'definition': 'A lipid similar in structure to a '
                                        'glycerolipid but in which one or more '
                                        'of the carbon atoms on glycerol is '
                                        'bonded to an alkyl chain via an ether '
                                        'linkage, as opposed to the usual '
                                        'ester linkage.',
                          'parents': ['CHEBI:18059', 'CHEBI:52575']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 10-11: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}