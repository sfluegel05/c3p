"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    This function checks for the presence of metal ions and common anions found in mineral nutrients.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define common metal ions in mineral nutrients
    metal_ions = {'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',  # Alkali metals
                  'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',  # Alkaline earth metals
                  'Al', 'Ga', 'In', 'Tl',              # Group 13 metals
                  'Zn', 'Cd', 'Hg',                    # Group 12 metals
                  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Ag', 'Au',  # Transition metals
                  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'}

    # Define common anions in mineral nutrients
    common_anions = {'Cl', 'Br', 'I',             # Halides
                     'O', 'S', 'P', 'N', 'Si'}    # Oxides, sulfates, phosphates, etc.
    
    # Check for presence of metal ions
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not any(metal in atoms for metal in metal_ions):
        return False, "No metal ion typically found in mineral nutrients"
    
    # Check for presence of common anions
    if not any(anion in atoms for anion in common_anions):
        return False, "No common anion typically found in mineral nutrients"
    
    return True, "Contains metal ion and anion typical of mineral nutrients"