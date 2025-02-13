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
    
    # Define SMARTS pattern for glycerol backbone with OH at 2,3-positions and 1-ester linkage
    # Simple ester linkage at position 1
    glycerol_ester_1_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)C")
    glycerol_ester_1_pattern_racemic = Chem.MolFromSmarts("OC[C@@H](O)COC(=O)C")
    
    if not mol.HasSubstructMatch(glycerol_ester_1_pattern) and not mol.HasSubstructMatch(glycerol_ester_1_pattern_racemic):
        return False, "No ester linkage found at 1-position of glycerol"

    # Ensure only one ester linkage is present
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O"))
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"
    
    # Check for the presence of a fatty acid chain (long hydrocarbon chain)
    # Consider a flexible definition for a long carbon chain
    long_chain_pattern = Chem.MolFromSmarts("C(=O)O[CX4][CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No fatty acid chain (long hydrocarbon chain) found"

    return True, "Contains glycerol backbone with one fatty acid chain attached via ester bond at position 1"