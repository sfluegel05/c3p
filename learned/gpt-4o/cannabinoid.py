"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    A cannabinoid is identified by characteristic structural motifs and functional groups, 
    common in pharmacologically active compounds from or resembling the Cannabis plant.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broaden pattern for cannabinoid-like core structures
    # Example: Patterns for aromatic rings with common attachments like oxygens 
    cannabinoid_core_patterns = [
        Chem.MolFromSmarts('c1cc(O)c(C)c(C=CC(=O)O)c1'),  # General THC-like structure
        Chem.MolFromSmarts('c1cc(O)c(CO)cc1'),           # CBD-like
        Chem.MolFromSmarts('C=Cc1cc(O)c(C)cc1'),         # Delta-9-like core
    ]

    # Check for any one of the cannabinoid core patterns
    for pattern in cannabinoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Cannabinoid-like core structure detected"
    
    # Additional common characteristic features (presence of an ether, ester, or amide group)
    functional_groups = [
        Chem.MolFromSmarts('CO'),                 # Ethers
        Chem.MolFromSmarts('OC(=O)'),             # Esters
        Chem.MolFromSmarts('NC(=O)'),             # Amides
    ]

    for pattern in functional_groups:
        if mol.HasSubstructMatch(pattern):
            return True, "Interesting functional group indicative of cannabinoids detected"
    
    return False, "No cannabinoid characteristic patterns detected"