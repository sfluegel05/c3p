"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride typically has a single fatty acid chain attached to the primary 
    hydroxyl group of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Glycerol backbone with ester group at primary position
    glycerol_ester_pattern = Chem.MolFromSmarts("O[C@@H]([C@@H](O)CO)C(=O)")
    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No ester linkage found at primary 1-position of glycerol"
    
    # Ensure ester linkage leads to a long carbon chain pattern (C3-C22)
    long_chain_pattern = Chem.MolFromSmarts("C(=O)OC[C](C)[O]~[CH2,CH](~[CH2,CH]){2,20}")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Chain too short, does not match expected fatty acid chain length"
    
    return True, "Contains glycerol backbone with an ester linkage at the 1-position and a long fatty acid chain"