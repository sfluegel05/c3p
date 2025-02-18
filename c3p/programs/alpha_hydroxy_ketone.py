"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:2468 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone has a ketone group (C=O) with a hydroxyl (-OH) on the adjacent (alpha) carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-hydroxy ketone: ketone (C=O) with adjacent OH-bearing carbon
    # Pattern matches [CX3](=O)-[C][OH] where:
    # - [CX3] is the ketone carbonyl carbon
    # - [C] is the alpha carbon
    # - [OH] is the hydroxyl group on the alpha carbon
    pattern = Chem.MolFromSmarts('[CX3](=O)-[C][OH]')
    
    if mol.HasSubstructMatch(pattern):
        return True, "Contains ketone with hydroxyl group on adjacent alpha-carbon"
    else:
        return False, "No ketone with alpha-hydroxy group detected"