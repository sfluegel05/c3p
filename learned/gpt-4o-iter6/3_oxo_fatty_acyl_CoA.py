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

    # Define a comprehensive pattern for the CoA moiety
    coa_core_pattern = Chem.MolFromSmarts("[C@@H]1N(C=NC2=C1N=CN=C2N)C3C(C(C(O3)(COP(=O)(O)OCC[N+](C)(C)C))O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_core_pattern):
        return False, "No core CoA moiety found"

    # Define a pattern for the 3-oxo group with flexibility
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    if not mol.HasSubstructMatch(oxo_fatty_acid_pattern):
        return False, "No 3-oxo-fatty acid group found"
    
    # Look for the thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Structure matches 3-oxo-fatty acyl-CoA with core CoA moiety, 3-oxo group, and thioester linkage"

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