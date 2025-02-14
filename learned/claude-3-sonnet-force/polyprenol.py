"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:27480 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has the general formula H-[CH2C(Me)=CHCH2]nOH, where the carbon skeleton
    is composed of more than one isoprene unit.

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
    
    # Look for isoprene unit pattern (C=C-C-C=C or C=C-C=C-C)
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH2][CH2][CH2]=[C]|[CH2]=[C]([CH3])[CH]=[CH][CH2]=[C]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Check for at least two isoprene units
    if len(isoprene_matches) < 2:
        return False, "Less than two isoprene units found"
    
    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_match = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_match:
        return False, "No hydroxyl group found"
    
    # Check for a hydrocarbon skeleton with possible additional oxygen atoms
    heteroatom_pattern = Chem.MolFromSmarts("[!#6;!#1;!#8]")  # Non-carbon, non-hydrogen, non-oxygen atoms
    heteroatom_matches = mol.GetSubstructMatches(heteroatom_pattern)
    if len(heteroatom_matches) > 0:
        return False, "Found heteroatoms other than oxygen"
    
    return True, "Molecule has the general formula H-[CH2C(Me)=CHCH2]nOH with a hydrocarbon skeleton"