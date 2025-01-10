"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.
    A para-terphenyl consists of a 1,4-diphenylbenzene skeleton, potentially with substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule resembles a para-terphenyl structure, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Modified SMARTS pattern for para-terphenyl
    para_terphenyl_pattern = Chem.MolFromSmarts("c1(c2ccc(cc2)c3ccccc3)ccc(cc1)-c4cccc5ccccc45")
    
    # Check if the molecule matches the para-terphenyl pattern
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains a 1,4-diphenylbenzene skeleton characteristic of para-terphenyl"
    else:
        return False, "Does not contain a distinct 1,4-diphenylbenzene skeleton"