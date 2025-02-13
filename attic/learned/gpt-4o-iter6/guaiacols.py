"""
Classifies: CHEBI:134251 guaiacols
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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phenol pattern: benzene ring with an -OH group
    phenol_pattern = Chem.MolFromSmarts('c1cc(O)ccc1')
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol group found"

    # Look for ortho-methoxy pattern: methoxy at ortho position to -OH group
    guaiacol_pattern = Chem.MolFromSmarts('c1cc(OC)c(O)cc1')
    if not mol.HasSubstructMatch(guaiacol_pattern):
        return False, "No ortho-methoxy group found"

    return True, "Molecule is classified as a guaiacol: contains phenol with ortho-position methoxy group"