"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: CHEBI:77518 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester has a butyryl group (CH2CH2CH2CO-) attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define butyrate ester pattern: four carbons (including carbonyl) followed by ester oxygen
    # SMARTS: CCCC(=O)O where oxygen is connected to any atom (ester bond)
    pattern = Chem.MolFromSmarts("CCCC(=O)[OX2]")
    
    # Check for presence of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains butyryl group (CCCC=O) connected via ester bond"
    else:
        return False, "No butyrate ester group detected"