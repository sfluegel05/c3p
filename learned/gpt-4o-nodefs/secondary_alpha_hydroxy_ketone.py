"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone is defined by a ketone group with an adjacent secondary alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised SMARTS pattern for secondary alpha-hydroxy ketone
    # C[C@H](O)[C](C)=O where [C@H] can be any stereochemistry due to '@' chiral marker
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CH2X4][C@H](O)[CD3H1](=O)")

    # Check for the secondary alpha-hydroxy ketone pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains secondary alpha-hydroxy ketone structure"

    # The pattern checks for a secondary (2-carbon) hydroxylated carbon next to a ketone.
    return False, "Does not contain secondary alpha-hydroxy ketone structure"