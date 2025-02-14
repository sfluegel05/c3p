"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime derived from an aliphatic aldehyde,
    meaning it contains an aldoxime group (R-CH=N-OH) and has no aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aromatic atoms
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
    
    # Define aldoxime SMARTS pattern
    # Carbon with one hydrogen, double-bonded to nitrogen, single-bonded to oxygen with hydrogen
    aldoxime_pattern = Chem.MolFromSmarts("[CH]=N[OH]")
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aldoxime functional group found"
    
    return True, "Contains aldoxime group and is aliphatic"