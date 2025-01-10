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

    # Refined Coenzyme A structure: including adenine ring, riboses, phosphates, and thiol linkage.
    coenzyme_a_pattern = Chem.MolFromSmarts("NC1=NC=CN=C1N2C=NC(N)N=C2.C3(O[C@H](COP(=O)(O)OP(=O)(O)OCC4O[C@@H](O)[C@H](O)[C@H]4O)O3).SC(=O)[C@@H](N)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)O")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Incomplete CoA structure not detected"
    
    # Thioester linkage pattern with some flexibility in the acyl-CoA linkage.
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage with CoA not found"
    
    # Fatty acyl chain pattern allowing variations in length and double bonds.
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]{4,}")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficient hydrocarbon chain detected"

    return True, "Contains complete Coenzyme A structure, thioester linkage, and fatty acyl moiety"