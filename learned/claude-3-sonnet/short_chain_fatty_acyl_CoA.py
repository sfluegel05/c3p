"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:35601 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("C1=NC2=NC(=NC(=N2)N)N1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OCC4OC(N5C=NC6=C5N=CN=C6N)C(O)C4O)O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Look for fatty acid moiety pattern (short carbon chain with acid group)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX4][CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 1:
        return False, f"Found {len(fatty_acid_matches)} fatty acid moieties, expected 1"
    
    # Check if fatty acid moiety is connected to CoA via thioester bond
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester bonds, expected 1"
    
    # Count carbon atoms in the fatty acid chain
    chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetIsInRingSize() == 0)
    if chain_length < 4 or chain_length > 6:
        return False, "Fatty acid chain length not in the range of 4-6 carbons"
    
    return True, "Contains a CoA moiety connected to a short-chain fatty acid via a thioester bond"