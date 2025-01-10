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
    
    # Define a broader range of metal and metalloid elements found in minerals
    known_metals = {
        'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Cs', 'Ba', 'La', 'Ce',
        'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Th', 'Pa', 'U'
    }
    
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
        Chem.MolFromSmarts('[Cl-]'),  # Chloride
        Chem.MolFromSmarts('B([O-])[O-]'),  # Borate
        Chem.MolFromSmarts('Si(=[O-])(O[Si]([O-])=O)[O-]')  # Silicate
    ]
    
    for pattern in anion_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains characteristic anions and metal elements typical of minerals"
    
    # Specific patterns for unique mineral identifications
    specifics = [
        Chem.MolFromSmarts('[Fe++].[S-][S-]'),  # Pyrite
        Chem.MolFromSmarts('[S--].[Fe+3].[As-]'),  # Arsenopyrite
        Chem.MolFromSmarts('[Ni]=S=[Ni]=S=[Ni]'),  # Heazlewoodite
        Chem.MolFromSmarts('Cl[O-].[Ca+2].Cl[O-]'),  # Calcium hypochlorite
        Chem.MolFromSmarts('O=[Si]([O-])O[Si](=O)[O-].[Al+3].[Al+3]')  # Kaolinite
    ]
    
    for specific in specifics:
        if mol.HasSubstructMatch(specific):
            return True, "Contains crystalline structure typical of some minerals"
    
    # Check for common hydrate structures based on logic
    if metal_count > 0:
        return True, "Metal elements suggest potential mineral classification"
    
    return False, "Does not meet typical mineral structure criteria"