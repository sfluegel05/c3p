"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the class of para-terphenyl compounds based on its SMILES string.
    A para-terphenyl is characterized by a structure consisting of three phenyl rings, potentially with 
    substitutions, where two outer phenyl groups are connected via a central phenyl group in a linear fashion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a para-terphenyl or derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible SMARTS pattern
    # This pattern allows for typical para-terphenyl core connectivity with potential substitutions
    para_terphenyl_pattern = Chem.MolFromSmarts('c1(cccc1)-c2ccccc2-c3ccccc3')
    
    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        return False, "Does not match para-terphenyl core structure"
    
    # Additional checks or descriptors can be added here to refine classification, 
    # such as expected functional groups or specific bonds consistent with known derivatives

    return True, "Contains para-terphenyl core structure"