"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA contains a fatty acyl chain with one carbon-carbon double bond and a CoA group.

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
    
    # Define CoA backbones substructure
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA group found"
    
    # Identify long carbon chain with one double bond
    # General pattern: carbon chains with double bonds (C=C)
    unsat_fatty_acid_pattern = Chem.MolFromSmarts("C=CC")
    double_bond_matches = mol.GetSubstructMatches(unsat_fatty_acid_pattern)
    
    if len(double_bond_matches) != 1:
        return False, f"Found {len(double_bond_matches)} double bonds, need exactly 1"
    
    return True, "Contains CoA and a fatty acyl chain with one carbon-carbon double bond"

# Example usage:
# smiles = "CCCCC/C=C\\CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
# result, reason = is_monounsaturated_fatty_acyl_CoA(smiles)
# print(result, reason)