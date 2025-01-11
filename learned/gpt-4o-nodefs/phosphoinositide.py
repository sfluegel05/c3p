"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide must contain a myo-inositol ring with phosphate groups and fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for myo-inositol, phosphate, and glycerol backbone
    inositol_pattern = Chem.MolFromSmarts("C1([C@@H]([C@H]([C@@H](C([C@@H]([C@@H]1O)O)O)O)O)O)")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    glycerol_pattern = Chem.MolFromSmarts("OC[C@@H](O)CO")
    
    # Check for myo-inositol structure
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol structure found"
    
    # Check for at least one phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"
    
    # Check for a glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for at least one fatty acid chain (2 or more carbons)
    # Assuming a simple representation here as long carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("CCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No fatty acid chains found"
    
    # If all checks pass
    return True, "Contains myo-inositol structure, phosphate group(s), glycerol backbone, and fatty acid chains"

# Example usage
# smiles_str = "a correct SMILES string for a known phosphoinositide"
# is_phosphoinositide(smiles_str)