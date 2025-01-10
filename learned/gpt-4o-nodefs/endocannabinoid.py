"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically have a long polyunsaturated alkyl chain and may
    include functionalities like ethanolamide or glycerol derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flexible pattern for polyunsaturated alkyl chain (at least one double bond)
    unsaturated_chain_pattern = Chem.MolFromSmarts("C=C.CCC.C=C")  # indicative of non-conjugated double bonds
    if not mol.HasSubstructMatch(unsaturated_chain_pattern):
        return False, "No long polyunsaturated alkyl chain found"

    # Check for either an ethanolamide group or a glycerol derivative
    ethanolamide_pattern = Chem.MolFromSmarts("NC(=O)CO")  # N-acyl ethanolamine
    glycerol_pattern = Chem.MolFromSmarts("OC(CO)CO")  # Simple glycerol structure
    if not (mol.HasSubstructMatch(ethanolamide_pattern) or mol.HasSubstructMatch(glycerol_pattern)):
        return False, "No recognizable group (ethanolamide or glycerol derivative) found"

    # Counts number of carbon atoms for potential long chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if c_count < 15:
        return False, "Too few carbons for typical endocannabinoid"

    return True, "Matches endocannabinoid characteristics with polyunsaturated alkyl chain and functional group"