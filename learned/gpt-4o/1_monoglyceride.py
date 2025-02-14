"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride contains a glycerol backbone with one fatty acid chain
    attached via an ester bond at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (OH groups on 1, 2, 3 carbons)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ester linkage at the 1-position
    # -OH replaced with an ester group at the start of the glycerol
    ester_pattern_1_position = Chem.MolFromSmarts("OC(=O)C(O)C(O)")
    if not mol.HasSubstructMatch(ester_pattern_1_position):
        return False, "No ester linkage found at 1-position of glycerol"
    
    # Ensure only one ester linkage is present
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O"))
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"
    
    # Check for the presence of fatty acid chains
    long_chain_pattern = Chem.MolFromSmarts("C(CCCCCCCC)C(=O)O")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long-chain characteristic of fatty acid"
    
    return True, "Contains glycerol backbone with one fatty acid chain attached via ester bond at position 1"