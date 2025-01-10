"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string.
    A mineral is generally characterized by specific cation-anion structures, often having metal elements and specific common anions.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a mineral, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define potential metal and metalloid elements and common anions in minerals
    metal_elements = {'Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Ba', 'Cs', 'La', 'Zn', 'Ni', 'Sb'}
    metalloids_oxides = {'B', 'Si', 'As', 'P', 'Al'}
    common_anions = {'Cl', 'F', 'O', 'S'}
    
    # Check for metal elements - required for mineral classification
    metal_count = sum(atom.GetSymbol() in metal_elements for atom in mol.GetAtoms())
    metalloid_oxide_count = sum(atom.GetSymbol() in metalloids_oxides for atom in mol.GetAtoms())
    anion_count = sum(atom.GetSymbol() in common_anions for atom in mol.GetAtoms())

    if metal_count < 1:
        return False, "No significant metal element typical of minerals detected"

    # Define anion patterns common in minerals
    anion_patterns = [
        Chem.MolFromSmarts('[O-]'),  # Oxide
        Chem.MolFromSmarts('P(=O)([O-])[O-]'),  # Phosphate
        Chem.MolFromSmarts('S(=O)([O-])[O-]'),  # Sulfate
        Chem.MolFromSmarts('[C](=O)([O-])'),  # Carbonate
        Chem.MolFromSmarts('[F-]'),  # Fluoride
        Chem.MolFromSmarts('[Cl-]'),  # Chloride
    ]
    
    for pattern in anion_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains characteristic anions and metal elements typical of minerals"

    # Check for hydrous minerals or complex mineral ions common in geology
    water_pattern = Chem.MolFromSmarts('O')
    water_matches = len(mol.GetSubstructMatches(water_pattern))
    if water_matches >= 3 and metal_count > 0:
        return True, "Contains multiple water molecules with metal elements, characteristic of hydrated minerals"
    
    # Specific checks for crystalline or unique mineral compounds that can be confirmed
    specific_mineral_patterns = [
        Chem.MolFromSmarts('[Fe++].[S-][S-]'),  # Pyrite pattern
        Chem.MolFromSmarts('[Ni]=S=[Ni]=S=[Ni]'),  # Heazlewoodite pattern
        Chem.MolFromSmarts('[S--].[Fe+3].[As-]')  # Arsenopyrite pattern
    ]
    for pattern in specific_mineral_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains crystalline structure typical of some minerals"
    
    return False, "Does not meet typical mineral structure criteria"