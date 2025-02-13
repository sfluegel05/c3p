"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    Furanocoumarins typically include a furochromone moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Furanocoumarin core pattern: furo[2,3-h]chromen or furo[3,2-g]chromen
    furanocoumarin_pattern = Chem.MolFromSmarts('O=c1cc2oc3ccccc3c2c(=O)[cH]c1')  # basic pattern for furanocoumarin

    if mol.HasSubstructMatch(furanocoumarin_pattern):
        return True, "Contains a core furochromone moiety typical of furanocoumarins"
    
    return False, "Does not contain key furanocoumarin structures"