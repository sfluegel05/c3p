"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a thioester linkage C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found between fatty acid and CoA"

    # Try to identify a long carbon chain (fatty acyl) pattern allowing flexibility
    # Acknowledging the possibility of unsaturation, branching, or other functional groups
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)[CX4;H2,CX3][CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No suitable fatty acyl chain found"

    # Define a more comprehensive pattern for coenzyme A
    coa_pattern = Chem.MolFromSmarts("NCC(=O)[C@H](O)C(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    return True, "Contains thioester linkage with fatty acyl chain and CoA structure"