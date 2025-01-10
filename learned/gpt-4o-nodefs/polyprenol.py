"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are oligomers with three or more isoprene units, often ending in
    an alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for general isoprene unit: CH2=C-CH2-CH=C
    isoprene_smarts = "[CH2]=[CH]-[CH2]-[CH]=[CH2]"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)

    # Find substructure matches
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    # Check for three or more isoprene units
    if len(isoprene_matches) < 3:
        return False, f"Only {len(isoprene_matches)} isoprene units found, at least 3 required"
    
    # Check for alcohol group (OH) at the terminal
    # We consider the possibility of various endings for polyprenols
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")  # Carbon with single-bonded OH
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Molecule does not end with an alcohol group"

    return True, "Molecule is a polyprenol with three or more isoprene units ending in an alcohol group"