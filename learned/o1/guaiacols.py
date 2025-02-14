"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    
    A guaiacol is defined as any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the guaiacol SMARTS pattern
    guaiacol_pattern = Chem.MolFromSmarts('c1c([OH])c([OCH3])cccc1')
    
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Molecule is a guaiacol (phenol with ortho methoxy group)"
    else:
        return False, "Molecule does not contain guaiacol substructure"