"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string.
    A mineral is generally characterized by specific cation-anion structures usually involving metal elements.
    
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
    
    # Define possible metal and metalloid elements found in minerals
    known_metals = {'Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Ba', 'Cs', 'La', 'Zn', 'Ni', 'Sb', 'Pd'}
    
    # Check if any metal is in the compound
    metal_count = sum(atom.GetSymbol() in known_metals for atom in mol.GetAtoms())
    if metal_count == 0:
        return False, "No significant metal element typical of minerals detected"
    
    # Define complex mineral anion patterns
    anion_patterns = [
        Chem.MolFromSmarts('[O-]'),  # Oxide
        Chem.MolFromSmarts('P(=O)([O-])[O-]'),  # Phosphate
        Chem.MolFromSmarts('S(=O)([O-])[O-]'),  # Sulfate
        Chem.MolFromSmarts('[C](=O)([O-])'),  # Carbonate
        Chem.MolFromSmarts('[F-]'),  # Fluoride
        Chem.MolFromSmarts('[Cl-]')  # Chloride
    ]
    
    for pattern in anion_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains characteristic anions and metal elements typical of minerals"
    
    # Check for common hydrate structures by counting oxygen atoms (to detect water presence)
    water_pattern = Chem.MolFromSmarts('O')
    water_count = len(mol.GetSubstructMatches(water_pattern))
    if water_count > 2 and metal_count > 0:
        return True, "Contains multiple water molecules with metal elements, characteristic of hydrated minerals"

    # Specific patterns for unique mineral identifications
    specifics = [
        Chem.MolFromSmarts('[Fe++].[S-][S-]'),  # Pyrite
        Chem.MolFromSmarts('[Ni]=S=[Ni]=S=[Ni]'),  # Heazlewoodite
        Chem.MolFromSmarts('[S--].[Fe+3].[As-]')  # Arsenopyrite
    ]
    
    for specific in specifics:
        if mol.HasSubstructMatch(specific):
            return True, "Contains crystalline structure typical of some minerals"
    
    return False, "Does not meet typical mineral structure criteria"