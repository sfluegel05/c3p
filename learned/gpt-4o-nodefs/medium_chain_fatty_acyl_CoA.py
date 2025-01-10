"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A moiety pattern
    coa_pattern = "[C@H](C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C(C(N)=O)=O)C(N(CCSC)=O)=O)O"
    coa_mol = Chem.MolFromSmarts(coa_pattern)
    if not mol.HasSubstructMatch(coa_mol):
        return False, "No CoA moiety found"

    # Thioester bond pattern
    thioester_pattern = "C(=O)SC"
    thioester_mol = Chem.MolFromSmarts(thioester_pattern)
    if not mol.HasSubstructMatch(thioester_mol):
        return False, "No thioester linkage found"

    # Identify the fatty acyl chain length
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    chain_lengths = [sum(1 for atom in match if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6) for match in carbon_chain_matches]
    
    # Check if it contains a medium chain (6-12 carbons)
    if any(6 <= length <= 12 for length in chain_lengths):
        return True, "Contains a CoA moiety with medium-chain fatty acyl group"
    
    return False, "Fatty acyl chain not of medium length"