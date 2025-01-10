"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is characterized by the presence of a hydroxyl group adjacent to a ketone group, defined as R1-CH(OH)-C(=O)-R2, where R1 and R2 can vary widely.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Update SMARTS pattern to catch more configurations of alpha-hydroxy ketones
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[C;R0,R1;H1,H2](O)[C](=O)")

    # Check if the molecule contains the alpha-hydroxy ketone substructure pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains an alpha-hydroxy ketone group"
    else:
        return False, "Does not contain an alpha-hydroxy ketone group"