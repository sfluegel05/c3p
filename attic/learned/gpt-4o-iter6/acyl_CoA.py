"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of a carboxylic acid and coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """  
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thioester bond (S-C(=O)-C)
    thioester_pattern = Chem.MolFromSmarts("SCC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for Coenzyme A backbone
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(OP(O)(O)=O)(OP(O)(O)=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    return True, "Contains thioester linkage and Coenzyme A structure"

# Example usage:
# smiles_example = "S(C(=O)CCC)CCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
# result, reason = is_acyl_CoA(smiles_example)
# print(result, reason)