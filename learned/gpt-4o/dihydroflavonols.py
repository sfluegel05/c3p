"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a generalized SMARTS pattern for dihydroflavonol
    # This encompasses the chromanone core with a hydroxyl group at C3
    dihydroflavonol_pattern = Chem.MolFromSmarts("C1C(O)OC2=CC=CC=C2C1=O")
    
    # Ensure a hydroxyl (O) is connected to the C3 carbon 
    hydroxyl_at_C3_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H]C(=O)Cc2ccccc12")
    
    # Check if the molecule matches either pattern
    if mol.HasSubstructMatch(dihydroflavonol_pattern) or mol.HasSubstructMatch(hydroxyl_at_C3_pattern):
        return True, "Contains dihydroflavonol structure"

    return False, "Does not fit dihydroflavonol pattern adequately"