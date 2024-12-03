"""
Classifies: CHEBI:76179 O-(S-acylpantetheine-4'-phosphoryl)serine(1-) residue
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_O__S_acylpantetheine_4__phosphoryl_serine_1___residue(smiles: str):
    """
    Determines if a molecule is an O-(S-acylpantetheine-4'-phosphoryl)serine(1-) residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-(S-acylpantetheine-4'-phosphoryl)serine(1-) residue, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of specific substructures
    substructures = [
        "C(C(COP([O-])(=O)OC[C@H](N*)C(=O)*)C)(C)O",  # Pantetheine substructure
        "C(=O)NCCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OC[C@H](N*)C(=O)*",  # Serine substructure
        "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OC[C@H](N*)C(=O)*"  # Acylpantetheine substructure
    ]

    for substructure in substructures:
        if not mol.HasSubstructMatch(Chem.MolFromSmarts(substructure)):
            return False, f"Missing substructure: {substructure}"

    # Check for the anionic state by looking for the phospho group with a negative charge
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("P([O-])(=O)OC")):
        return False, "Missing anionic phospho group"

    return True, "Molecule is an O-(S-acylpantetheine-4'-phosphoryl)serine(1-) residue"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76179',
                          'name': "O-(S-acylpantetheine-4'-phosphoryl)serine(1-) "
                                  'residue',
                          'definition': 'An anionic amino-acid residue formed '
                                        'by proton loss from the phospho group '
                                        'of an '
                                        "O-(S-acylpantetheine-4'-phosphoryl)serine "
                                        'residue, the principal form present '
                                        'at pH 7.3.',
                          'parents': ['CHEBI:64898']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 13-14: malformed \\N character escape (<string>, line '
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