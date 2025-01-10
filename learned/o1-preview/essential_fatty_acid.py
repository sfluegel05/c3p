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
    An essential fatty acid is a polyunsaturated fatty acid with an absolute dietary requirement.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdchem
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)O)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No terminal carboxylic acid group found"
    
    # Check for long hydrocarbon chain
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 16:
        return False, f"Chain too short ({num_carbons} carbons), needs at least 16 carbons"
    
    # Count number of double bonds (C=C)
    double_bond = Chem.MolFromSmarts("C=C")
    num_double_bonds = len(mol.GetSubstructMatches(double_bond))
    if num_double_bonds < 2:
        return False, f"Not polyunsaturated, only {num_double_bonds} double bond(s) found"
    
    # Check if double bonds are cis
    cis_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo == rdchem.BondStereo.STEREOZ:
                cis_double_bonds += 1
    if cis_double_bonds < 2:
        return False, f"Less than 2 cis double bonds ({cis_double_bonds} found)"
    
    # Check if molecule is a fatty acid (no other functional groups)
    allowed_atoms = {1,6,8}  # H, C, O atoms only
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Atom {atom.GetSymbol()} not allowed in fatty acids"
        if atom.GetAtomicNum() == 8 and atom.GetDegree() > 1 and not atom.IsInRing():
            # Allow only carboxylic acid oxygens
            connected_atoms = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            if 6 not in connected_atoms:
                return False, "Oxygen atom not part of carboxylic acid group"
    
    return True, "Molecule is a polyunsaturated fatty acid with an absolute dietary requirement"

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
    'attempt': 2,
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