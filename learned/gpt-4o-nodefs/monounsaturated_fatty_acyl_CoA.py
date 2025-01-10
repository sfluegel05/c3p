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

    # Define a flexible pattern for an aliphatic chain with one double bond
    # The pattern '[C;H2]C=C[C;H2]' corresponds to any chain segment with a double bond and 2 hydrogens (indicative of the chain's ends)
    double_bond_pattern = Chem.MolFromSmarts('[C]!@!=[C]')
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)

    # Filter matches to count only those in long carbon chains
    valid_double_bond_count = 0
    for match in double_bond_matches:
        # Simple conditions to ensure it's part of a long carbon chain
        chain_length = len(Chem.GetMolFrags(Chem.PathToSubmol(mol, match, useQueryQueryMatches=True), asMols=False))
        if chain_length >= 8:  # Arbitrarily chosen, adjust if needed
            valid_double_bond_count += 1

    if valid_double_bond_count != 1:
        return False, f"Found {valid_double_bond_count} aliphatic double bonds, need exactly 1 for monounsaturation"

    return True, "Monounsaturated fatty acyl-CoA compound identified"

# Example metadata
__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Coenzyme A linked with a monounsaturated fatty acid chain.'
    }
}