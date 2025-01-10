"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA contains a fatty acyl chain with one carbon-carbon double bond 
    and a CoA group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA group substructure - includes phosphate and ribose moiety components
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA group found"

    # Define a pattern for carbon-carbon double bonds in aliphatic chains
    double_bond_pattern = Chem.MolFromSmarts("C=C")  # Flexible C=C bond for aliphatic chains
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)

    # Ensure there's exactly one carbon-carbon double bond in the fatty acyl chain
    num_double_bonds = len(double_bond_matches)
    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds in aliphatic chain, need exactly 1"

    # Check for the fatty acyl chain pattern starting from the CoA
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)[!#1]C=C")  # C=O followed by any carbon chain with a double bond
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
        return False, "No valid monounsaturated fatty acyl chain found"
    
    return True, "Contains CoA and a fatty acyl chain with one carbon-carbon double bond"