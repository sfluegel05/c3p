"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic compound essential for the human body,
    typically consisting of an inorganic anion and a metal cation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True, reason if classified as a mineral nutrient;
                    False, reason otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for metal cations in essential mineral nutrients
    metal_cations_smarts = "[K+],[Na+],[Ca+2],[Mg+2],[Fe+3],[Zn+2],[Al+3],[Ba+2],[Cs+],[La+3],[Sb+5]"
    metal_cations_pattern = Chem.MolFromSmarts(metal_cations_smarts)
    
    # Define SMARTS patterns for common anionic groups (sulfates, chlorides, etc.)
    anions_smarts = "[Cl-],[F-],[P]([O-])([O-])([O-])=O,[S]([O-])([O-])=O,[C]([O-])([O-])=O,[N]([O-])([O-])=O,[Si]([O-])([O-])([O-])"
    anions_pattern = Chem.MolFromSmarts(anions_smarts)
    
    # Check for presence of at least one metal cation
    if not mol.HasSubstructMatch(metal_cations_pattern):
        return False, "No essential metal cations found"
    
    # Check for presence of at least one anion or inorganic group
    if not mol.HasSubstructMatch(anions_pattern):
        return False, "No recognizable anionic groups found"
    
    # If both components are present, classify as a potential mineral nutrient
    return True, "Contains essential metal cation and anionic counterpart typical of mineral nutrients"