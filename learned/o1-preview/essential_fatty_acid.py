"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acid
Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a long-chain polyunsaturated fatty acid that cannot be synthesized by the body and must be obtained from the diet.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Identify the carboxyl carbon(s) to exclude
    carboxyl_carbons = set()
    for match in mol.GetSubstructMatches(carboxylic_acid):
        carboxyl_carbons.add(match[0])  # match[0] is the carbonyl carbon atom

    # Count total number of carbon atoms (excluding carboxyl carbon(s))
    total_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetIdx() not in carboxyl_carbons:
                total_carbons += 1

    if total_carbons < 16:
        return False, f"Chain length {total_carbons} is too short for essential fatty acid (minimum is 16 carbons)"

    # Count number of carbon-carbon double bonds and check stereochemistry
    double_bonds = 0
    cis_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                double_bonds += 1
                # Check stereochemistry
                stereo = bond.GetStereo()
                if stereo == Chem.rdchem.BondStereo.STEREOZ or stereo == Chem.rdchem.BondStereo.STEREOCIS:
                    cis_double_bonds += 1

    if double_bonds < 2:
        return False, f"Only {double_bonds} double bond(s) found, not enough for essential fatty acid (minimum is 2)"

    # Essential fatty acids typically have cis-configured double bonds
    # If stereochemistry is unspecified (STEREONONE), we cannot confirm cis configuration
    if cis_double_bonds < double_bonds:
        return False, f"Not all double bonds are in cis configuration (found {cis_double_bonds} cis out of {double_bonds})"

    # If all checks pass, classify as essential fatty acid
    return True, "Molecule is a long-chain polyunsaturated fatty acid with multiple cis double bonds"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'essential fatty acid',
        'definition': 'Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.',
        'parents': ['polyunsaturated fatty acid', 'fatty acid'],
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