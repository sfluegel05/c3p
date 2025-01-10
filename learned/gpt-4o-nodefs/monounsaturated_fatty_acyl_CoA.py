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
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA structure
    coa_pattern = Chem.MolFromSmarts('COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structure not found"

    # Check for exactly one double bond in the fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 1:
        return False, f"Found {len(double_bond_matches)} double bonds, need exactly 1 for monounsaturation"

    # Check for general aliphatic character of long-chain fatty acid
    # This is a hard problem, a simple check for length is a proxy
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # typically expect longer chains
        return False, "Not enough carbons for a fatty acid chain characteristic"

    return True, "Monounsaturated fatty acyl-CoA compound identified"

# Example metadata
__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Coenzyme A linked with a monounsaturated fatty acid chain.'
    }
}