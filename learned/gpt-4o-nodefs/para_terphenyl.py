"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the class of para-terphenyl compounds based on its SMILES string.
    A para-terphenyl is characterized by a structure consisting of three phenyl rings, where two outer phenyl 
    groups are connected via a central phenyl group in a linear fashion (para orientation).

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
    
    # Define an enhanced SMARTS pattern for para-terphenyl
    # Allow for potential substitutions at the para positions of each phenyl ring
    para_terphenyl_pattern = Chem.MolFromSmarts('c1(cc([#1,#6,#8,F,Cl,Br,S](c)*)cc1)-c2cc([#1,#6,#8,F,Cl,Br,S])ccc2-c3ccc([#1,#6,#8,F,Cl,Br,S])cc3')
    
    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        return False, "Does not match para-terphenyl core structure"
    
    # Check if substitutions adhere to known patterns in derivatives
    # This part can be extended to check for particular functional groups or combinations
    # that are characteristic of known para-terphenyl derivatives

    return True, "Contains para-terphenyl core structure"