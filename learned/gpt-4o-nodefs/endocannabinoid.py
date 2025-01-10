"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically have a long polyunsaturated alkyl chain and may
    contain an ethanolamide group.
    
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
    
    # Look for long alkyl chain with at least multiple double bonds (up to 4)
    polyunsaturated_pattern = Chem.MolFromSmarts("C=CCCC=CCCC=CCCC=CC")
    if not mol.HasSubstructMatch(polyunsaturated_pattern):
        return False, "No long polyunsaturated alkyl chain found"
    
    # Look for ethanolamide group
    ethanolamide_pattern = Chem.MolFromSmarts("NCCO")
    if not mol.HasSubstructMatch(ethanolamide_pattern):
        return False, "No ethanolamide group found"
    
    # Counts number of carbon atoms for long alkyl chain if needed
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if c_count < 20:
        return False, "Too few carbons for typical endocannabinoid"
    
    return True, "Matches endocannabinoid characteristics with polyunsaturated alkyl chain and ethanolamide group"