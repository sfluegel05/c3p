"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: CHEBI:27027 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is typically an inorganic compound containing metal cations and inorganic anions.

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

    # Check for the presence of metal cations (e.g., Ca, Mg, K, Na, etc.)
    metal_cations = {'Ca', 'Mg', 'K', 'Na', 'Ba', 'Cs', 'Zn', 'Fe', 'Al', 'La', 'Sb', 'Pd'}
    has_metal_cation = any(atom.GetSymbol() in metal_cations for atom in mol.GetAtoms())
    if not has_metal_cation:
        return False, "No metal cation found"

    # Check for the presence of inorganic anions (e.g., phosphate, sulfate, nitrate, etc.)
    inorganic_anions = {'O', 'P', 'S', 'N', 'Cl', 'F', 'Si'}
    has_inorganic_anion = any(atom.GetSymbol() in inorganic_anions for atom in mol.GetAtoms())
    if not has_inorganic_anion:
        return False, "No inorganic anion found"

    # Check if the molecule is inorganic (no carbon atoms or only carbon in carbonate-like structures)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) > 0:
        # Allow carbon only in carbonate-like structures (e.g., CO3)
        carbonate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])([OX2])")
        if not mol.HasSubstructMatch(carbonate_pattern):
            return False, "Contains carbon atoms not in carbonate-like structures"

    # Check if the molecule is a salt (contains both cations and anions)
    # This is a simple check based on the presence of charged atoms
    has_cation = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    has_anion = any(atom.GetFormalCharge() < 0 for atom in mol.GetAtoms())
    if not (has_cation and has_anion):
        return False, "Not a salt (missing cations or anions)"

    return True, "Contains metal cations and inorganic anions, typical of mineral nutrients"