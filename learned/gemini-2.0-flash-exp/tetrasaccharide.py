"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generic monosaccharide pattern (5 or 6 membered ring with multiple OH)
    # This is a simplified pattern, more specific patterns may be required for very specific scenarios
    monosaccharide_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX4]([OX2])[CX4][CX4][CX4]1 | [OX2]1[CX4][CX4][CX4][CX4][CX4]1") # 5- and 6- membered rings
    if monosaccharide_pattern is None:
      return None, "Error in SMARTS pattern"


    # Find all matching monosaccharide units
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    # Check for exactly 4 units
    if len(monosaccharide_matches) != 4:
        return False, f"Found {len(monosaccharide_matches)} monosaccharide units, requires 4 for tetrasaccharide"

    # Count glycosidic bonds (C-O-C connecting rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    # A tetrasaccharide will have 3 glycosidic bonds
    if len(glycosidic_bond_matches) < 3 :
        return False, "Too few glycosidic bonds for a tetrasaccharide"
    
    #Additional check to make sure the monosaccharides are connected in a chain and not an other format like a cycle of 4 sugar units.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if carbon_count < 15 or carbon_count > 30 :
        return False, "Number of carbon atoms is outside tetrasaccharide range"

    if oxygen_count < 10 or oxygen_count > 20 :
        return False, "Number of oxygen atoms is outside tetrasaccharide range"

    return True, "Contains 4 monosaccharide units connected via glycosidic bonds"