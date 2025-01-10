"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the class of para-terphenyl compounds based on its SMILES string.
    A para-terphenyl is characterized by a structure consisting of three phenyl rings in a linear connection.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for para-terphenyl
    # Matches three benzene rings connected linearly (para)
    para_terphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2-c3ccccc3')
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains para-terphenyl core structure"
    
    return False, "Does not match para-terphenyl core structure"