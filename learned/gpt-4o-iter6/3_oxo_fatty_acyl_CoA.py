"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    This involves checking for a CoA moiety, 3-oxo group, and a thioester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended pattern for the CoA moiety
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)COP(=O)(OCC1OC(CO1)n2cnc3c(ncnc23)N)CCNC(=O)C(C)C(=O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No core CoA moiety found"
    
    # More specific pattern for the 3-oxo group
    oxo_group_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo-fatty acid group found"

    # Checking for the thioester linkage 'C(=O)S'
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