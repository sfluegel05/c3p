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

    # Define a pattern for the 3-oxo-fatty acid group (ketone at position 3)
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CC")
    if not mol.HasSubstructMatch(oxo_fatty_acid_pattern):
        return False, "No 3-oxo-fatty acid group found"
    
    # Define a more flexible pattern for the core CoA moiety
    coa_core_pattern = Chem.MolFromSmarts("P(=O)(O)OC[C@@H]1O[C@H]([C@@H](OP(=O)(O)O)[C@H]1O)N2C=NC3=C2N=CN=C3N")
    if not mol.HasSubstructMatch(coa_core_pattern):
        return False, "No core CoA moiety found"

    # Look for the thioester linkage (C(=O)S-C, ester linkage with sulfur)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found between 3-oxo and CoA"

    return True, "Structure matches 3-oxo-fatty acyl-CoA with 3-oxo group, core CoA moiety, and thioester linkage"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:xxxx',  # Replace with the correct CHEBI ID for 3-oxo-fatty acyl-CoA
        'name': '3-oxo-fatty acyl-CoA',
        'definition': 'An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.'
    },
    'config': {
        # Configuration details for model testing, thresholds etc.
    },
    'message': None,
    'success': True,
    'error': '',
    'stdout': None,
}