"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:35681 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA has a fatty acid chain of 6 carbons or fewer attached to a coenzyme A (CoA) moiety via a thioester bond.

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
    
    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(C)(CO)C(C(=O)NCCC(=O)NCCS)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A (CoA) moiety found"
    
    # Look for thioester bond connecting fatty acid chain
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])SCC")
    thioester_match = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_match:
        return False, "No thioester bond found to connect fatty acid chain"
    
    # Find fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    chain_lengths = [len(match) for match in fatty_acid_matches]
    if not chain_lengths:
        return False, "No fatty acid chain found"
    if max(chain_lengths) > 6:
        return False, "Fatty acid chain is too long (>6 carbons)"
    
    # Check for disallowed functional groups on fatty acid chain
    disallowed_pattern = Chem.MolFromSmarts("[SX2,PX3,OX1]~[CX4,CX3]")
    disallowed_match = mol.GetSubstructMatches(disallowed_pattern)
    if disallowed_match:
        return False, "Fatty acid chain contains disallowed functional groups (S, P, O=O)"
    
    # Check for cyclic fatty acid chains
    cyclic_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    cyclic_match = mol.GetSubstructMatches(cyclic_pattern)
    if cyclic_match:
        return False, "Fatty acid chain is cyclic"
    
    return True, "Contains a fatty acid chain of 6 carbons or fewer attached to coenzyme A via a thioester bond"