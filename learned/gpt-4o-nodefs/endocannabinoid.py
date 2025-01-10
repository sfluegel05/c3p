"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids often have long polyunsaturated alkyl chains with functionalities 
    like ethanolamide or glycerol derivatives, excluding phosphate and large hydrophilic groupings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for polyunsaturated alkyl chains (non-conjugated double bonds)
    unsaturated_chain_pattern = Chem.MolFromSmarts("CCCC=CCCC=CCCC")
    if not mol.HasSubstructMatch(unsaturated_chain_pattern):
        return False, "No long polyunsaturated alkyl chain found"

    # Pattern for ethanolamide group
    ethanolamide_pattern = Chem.MolFromSmarts("N(CCO)C(=O)")
    # Pattern for glycerol derivative (common in endocannabinoids)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")

    # Check for recognizable functional groups
    if not (mol.HasSubstructMatch(ethanolamide_pattern) or 
            mol.HasSubstructMatch(glycerol_pattern)):
        return False, "No recognizable group (ethanolamide or glycerol derivative) found"

    # Exclude compounds with phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, likely not an endocannabinoid"

    # Ensure a reasonable number of carbon atoms for typical endocannabinoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbons for typical endocannabinoid"

    # If holds all checks, likely endocannabinoid
    return True, "Matches endocannabinoid characteristics with polyunsaturated alkyl chain and functional group"