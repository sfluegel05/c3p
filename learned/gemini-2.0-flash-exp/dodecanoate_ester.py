"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is any ester in which the carboxylic acid component is lauric acid (dodecanoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the dodecanoyl ester group.
    # Matches: CCCCCCCCCCCC(=O)O-[anything]
    dodecanoyl_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O[#6]") #The #6 means a carbon of any kind
    
    # Check if the pattern is present
    if not mol.HasSubstructMatch(dodecanoyl_ester_pattern):
       return False, "No dodecanoyl ester group detected"
    
    return True, "Contains at least one dodecanoate ester group."