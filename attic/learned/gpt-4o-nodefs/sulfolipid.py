"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonic/sulfate group pattern (OS(O)(=O)=O)
    sulfate_pattern = Chem.MolFromSmarts("[O-]S(=O)(=O)[O-]") 
    sulfate_matches = mol.HasSubstructMatch(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate group found"

    # Look for sugar-like structure by 5-6 membered ring with oxygens
    sugar_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)")
    sugar_matches = mol.HasSubstructMatch(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar-like structure found"
        
    # Look for long chain fatty acid structures
    long_chain_pattern = Chem.MolFromSmarts("C[CH2]{10,}")
    long_chain_matches = mol.HasSubstructMatch(long_chain_pattern)
    if not long_chain_matches:
        return False, "No long alkyl chain found"

    # If all key features are found
    return True, "Contains a sulfate group, sugar-like structure, and long alkyl chain"