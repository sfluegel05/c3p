"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Classifies if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA pattern recognizing key fragments of CoA
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)C[OP](=O)([O])O[O][C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"

    # Define thioester bond pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester bond not found"

    # Check for branched alkyl group pattern, e.g., aliphatic branch
    branch_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branched alkyl chain found"

    # Additional check could be added to verify chain length if needed

    return True, "Contains branched-chain fatty acyl-CoA structure"

# Example usage
smiles_example = "C[C@H](O)[C@H](C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
result, reason = is_branched_chain_fatty_acyl_CoA(smiles_example)
print(f"Is branched-chain fatty acyl-CoA: {result}, Reason: {reason}")