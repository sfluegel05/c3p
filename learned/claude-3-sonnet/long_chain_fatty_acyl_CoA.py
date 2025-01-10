"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    Long-chain fatty acyl-CoAs have:
    - CoA moiety
    - Thioester linkage
    - Fatty acid chain length C13-C22
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for adenine base (part of CoA)
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine moiety found (required for CoA)"

    # Check for ribose with phosphate (part of CoA) 
    ribose_phosphate = Chem.MolFromSmarts("OCC1OC(n2cnc3c2)C(O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "No ribose-phosphate moiety found (required for CoA)"

    # Check for pantothenic acid part (part of CoA)
    pantothenic = Chem.MolFromSmarts("NCCC(=O)NCCC(O)C(C)(C)")
    if not mol.HasSubstructMatch(pantothenic):
        return False, "No pantothenic acid moiety found (required for CoA)"

    # Check for thioester linkage
    thioester = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester linkage found"

    # Count carbons in the fatty acid chain
    # First, find the thioester carbon
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "Could not analyze fatty acid chain"
    
    # Count carbons in the main chain (excluding CoA)
    # This is an approximation - we count all carbons not part of CoA structure
    coa_carbons = 23  # Approximate number of carbons in CoA
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    chain_carbons = total_carbons - coa_carbons
    
    if chain_carbons < 13:
        return False, f"Fatty acid chain too short (C{chain_carbons}, need C13-C22)"
    if chain_carbons > 22:
        return False, f"Fatty acid chain too long (C{chain_carbons}, need C13-C22)"

    # Additional check for diphosphate bridge
    diphosphate = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)")
    if not mol.HasSubstructMatch(diphosphate):
        return False, "No diphosphate bridge found (required for CoA)"

    return True, f"Long-chain fatty acyl-CoA with approximately C{chain_carbons} fatty acid chain"