"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids often have long polyunsaturated or hydroxyl chain with 
    functionalities like ethanolamide or glycerol derivatives.

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
    
    # Define a more inclusive pattern for long polyunsaturated chains
    polyunsaturated_chain_pattern = Chem.MolFromSmarts("CCCC=CCCC=CCCC=CCCC")
    if not mol.HasSubstructMatch(polyunsaturated_chain_pattern):
        return False, "No long polyunsaturated alkyl chain found"

    # Patterns for ethanolamide and glycerol-derived groups
    ethanolamide_pattern = Chem.MolFromSmarts("NCC=O")
    glycerol_ether_pattern = Chem.MolFromSmarts("OCC(O)CO")

    # Check if at least one of these endocannabinoid functional groups is present
    functional_group_present = (
        mol.HasSubstructMatch(ethanolamide_pattern) or 
        mol.HasSubstructMatch(glycerol_ether_pattern)
    )
    if not functional_group_present:
        return False, "No recognizable endocannabinoid functional group found"

    # Filter out molecules with phosphate groups, relevant for misclassified entries
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, not typical for endocannabinoids"

    # Carbon count check for complexity ensuring it is at least 15 carbons, typical in endocannabinoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbons for a typical endocannabinoid"

    return True, "Matches characteristic endocannabinoid features with polyunsaturated chain and functional group"