"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide typically contains a myo-inositol ring with one or more phosphate groups
    and fatty acid chains.

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
    
    # Define SMARTS patterns for myo-inositol (consider hydroxy positioning and chiral centers)
    inositol_pattern = Chem.MolFromSmarts("C([C@H]1[C@H]([C@@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol structure found"

    # Check for phosphate groups (tripleS patterns to match P(=O)(O)(O) on inositols)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "Less than one phosphate group found"
    
    # Define pattern for glycerol backbone linking to inositol
    glycerol_with_fatty_acid = Chem.MolFromSmarts("OC[C@H](O)COC(=O)CC")
    if not mol.HasSubstructMatch(glycerol_with_fatty_acid):
        return False, "No glycerol backbone with fatty acid chains found"
    
    # Check for long hydrocarbon chains attached to glycerol backbone
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if len(mol.GetSubstructMatches(carbon_chain_pattern)) < 2:  # At least two long chains
        return False, "Not enough hydrocarbon chains found to match fatty acid chains"

    return True, "Contains myo-inositol ring with phosphate group(s) and fatty acid chains"

# Example usage
# smiles_str = "a SMILES string for a known phosphoinositide"
# is_phosphoinositide(smiles_str)