"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester — Any carboxylic ester where the carboxylic acid component is acetic acid.
This program checks whether the molecule contains an ester fragment of the type O-C(=O)-CH3.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is defined as any carboxylic ester where the carboxylic acid component is acetic acid,
    meaning the molecule must contain an ester group with an acyl part of CH3-C(=O)-.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an acetate ester, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the acetate ester moiety: O-C(=O)-CH3 
    acetate_pattern = Chem.MolFromSmarts("O-C(=O)[CH3]")
    
    # Check if the molecule has the acetate ester substructure
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Contains the acetate ester moiety (O-C(=O)-CH3)"
    else:
        return False, "Does not contain the acetate ester moiety (O-C(=O)-CH3)"