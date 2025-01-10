"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broad pattern for the CoA moiety, considering key groups
    coa_core_pattern = Chem.MolFromSmarts("COP(=O)(O)OC[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O")  # Simplified
    if not mol.HasSubstructMatch(coa_core_pattern):
        return False, "No core CoA moiety found"
    
    # Pattern for the 3-oxo-fatty acid group
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    if not mol.HasSubstructMatch(oxo_fatty_acid_pattern):
        return False, "No 3-oxo-fatty acid group found"

    # Look for the thioester linkage 'C(=O)S'
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Structure matches 3-oxo-fatty acyl-CoA with core CoA moiety, 3-oxo group, and thioester linkage"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25020',  # Hypothetical ID for 3-oxo-fatty acyl-CoA
        'name': '3-oxo-fatty acyl-CoA',
        'definition': 'An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.'
    },
    'config': {
        # Additional configuration information
    },
    'message': None,
    'success': True,
    'error': '',
    'stdout': None,
}