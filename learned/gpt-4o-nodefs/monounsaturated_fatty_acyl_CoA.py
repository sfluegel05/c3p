"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Convert SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pattern for Coenzyme A (CoA)
    coa_pattern = Chem.MolFromSmarts('COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structure not found"

    # Define pattern for one double bond in an aliphatic chain
    # An aliphatic chain here excludes functional groups like esters, amides, etc.
    # Here we are assuming a simple double bond with some assumptions about surrounding chain
    long_aliphatic_with_one_double_bond = Chem.MolFromSmarts('CCC=CCCC')  # Assuming a simple chain pattern

    # Check for exact count of double bonds within aliphatic chain context
    double_bond_matches = mol.GetSubstructMatches(long_aliphatic_with_one_double_bond)
    if len(double_bond_matches) != 1:
        return False, f"Found {len(double_bond_matches)} aliphatic double bonds, need exactly 1 for monounsaturation"

    return True, "Monounsaturated fatty acyl-CoA compound identified"

# Example metadata
__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Coenzyme A linked with a monounsaturated fatty acid chain.'
    }
}