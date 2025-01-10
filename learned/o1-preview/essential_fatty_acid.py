"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acid
Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Get the carbon chain (exclude functional groups)
    # Assume the longest carbon chain represents the fatty acid chain
    chains = mol.GetSubstructMatches(Chem.MolFromSmarts("[C;X4;!$(C(-O)(=O))]*"))
    if not chains:
        return False, "No carbon chain found"
    chain_length = max(len(chain) for chain in chains)

    # Essential fatty acids typically have chain lengths >= 16
    if chain_length < 16:
        return False, f"Chain length {chain_length} is too short for essential fatty acid"

    # Count the number of double bonds in the chain
    db_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(db_pattern))

    # Essential fatty acids are polyunsaturated with multiple double bonds
    if double_bonds < 2:
        return False, f"Only {double_bonds} double bond(s) found, not enough for essential fatty acid"

    # Check for cis configuration of double bonds
    # Cis double bonds have bond stereochemistry of type E or Z
    cis_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo == Chem.rdchem.BondStereo.STEREOZ:
                cis_bonds += 1

    if cis_bonds < double_bonds:
        return False, f"Not all double bonds are in cis configuration (found {cis_bonds} cis out of {double_bonds})"

    # If all checks pass, classify as essential fatty acid
    return True, "Molecule is a polyunsaturated fatty acid with long chain and multiple cis double bonds"

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