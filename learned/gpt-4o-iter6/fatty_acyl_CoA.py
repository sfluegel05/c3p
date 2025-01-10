"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA contains a Coenzyme A conjugated to a fatty acyl chain via a thioester linkage.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A structure including adenine ring, riboses, phosphates, and thiol.
    coenzyme_a_pattern = Chem.MolFromSmarts("N=C1N=CN(C2=C1N=CN2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OCC4C(C(C(O4)COP(=O)(O)O)O)O)O)O")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Incomplete Coenzyme A structure not detected"
    
    # Thioester pattern focusing on the acyl-CoA linkage.
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage with CoA not found"
    
    # Fatty acyl chain pattern allowing variations in length and double bonds.
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]=[CX4,CX3][CX4,CX3]~[CX4,CX3]{8,}")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficient hydrocarbon chain detected"

    return True, "Contains complete Coenzyme A structure, thioester linkage, and fatty acyl moiety"