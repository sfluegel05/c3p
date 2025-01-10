"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has a single fatty acid chain attached to the primary hydroxyl group of glycerol.

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

    # Look for glycerol backbone pattern with two free hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two free hydroxyl groups found"

    # Pattern matching for ester linkage at the primary oxygen of glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)OCC(O)CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Ester linkage not found at primary oxygen of glycerol"

    # Ensure the ester linkage leads to a long carbon chain (at least 3 carbon atoms)
    # Start after the ester carbonyl [C(=O)] leading to a long alkyl chain
    long_chain_pattern = Chem.MolFromSmarts("C(=O)O[C]"+"~[CH2~CH/*]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Ester linkage does not lead to a sufficiently long carbon chain"

    return True, "Contains glycerol backbone with a single fatty acid chain esterified at the 1-position"