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

    # Define the CoA group substructure (including possible variations in representations)
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CC(N)C(=O)S")

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA group found"

    # Define a pattern for carbon-carbon double bonds in the main chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    
    # Get matches for double bonds
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)

    # Further ensure that the double bond(s) belong to a true fatty acyl chain
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C")
    
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
        return False, "No valid fatty acyl chain found"

    # Ensure there's exactly one carbon-carbon double bond within that fatty acyl chain
    num_double_bonds = len(double_bond_matches)
    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds in main chain, need exactly 1"

    return True, "Contains CoA and a fatty acyl chain with one carbon-carbon double bond"

# Example usage:
# smiles = "CCCCC/C=C\\CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
# result, reason = is_monounsaturated_fatty_acyl_CoA(smiles)
# print(result, reason)