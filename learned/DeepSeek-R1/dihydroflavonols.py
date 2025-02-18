"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:???? dihydroflavonols (3-hydroxyflavanones)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for 3-hydroxyflavanone core
    # The pattern matches the flavanone skeleton with a hydroxyl group on position 3
    pattern = Chem.MolFromSmarts("[OH]C1C(=O)C2=CC=CC=C2C1O")
    
    if not mol.HasSubstructMatch(pattern):
        return False, "Does not contain the 3-hydroxyflavanone core"
    
    return True, "Contains 3-hydroxyflavanone structure"