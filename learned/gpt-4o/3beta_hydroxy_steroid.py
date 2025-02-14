"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the basic steroid backbone pattern (ABCD ring system)
    # R - cyclic carbon pattern, Cn - junction carbon with stereocenters flexibility (using @? for unspecified stereochemistry)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3(C)CCC4')
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Define the 3beta-hydroxy pattern with beta orientation at C3
    # @? indicates we are looking for substituents at stereocenters without strict stereochemical constraints
    hydroxy_3beta_pattern = Chem.MolFromSmarts('[C@]1([*:2])C2=C(C[C@]3([C@@]2CC[C@@H]3[*:3])[C@@](C1)([*:1])[#6])O')
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "3beta-hydroxy group not found"
    
    return True, "Contains steroid backbone and a 3beta-hydroxy group"