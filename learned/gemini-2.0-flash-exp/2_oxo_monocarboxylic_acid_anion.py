"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion has a C=O group at the 2-position relative to a carboxylic acid anion (C(=O)[O-]) group.
    It has exactly one carboxylic acid group

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Canonicalize the molecule
    smiles = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles)
    
    # Define carboyxlic acid pattern (allowing both -O- and -OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX1H0-]") #C(=O)[O-] or C(=O)OH
    
    # Define the 2-oxo-monocarboxylic acid pattern
    oxo_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][CX4](=O)[OX1H0-]")

    # Check if molecule has exactly one carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Molecule contains {len(carboxylic_acid_matches)} carboxylic acid groups, expected exactly 1"
        
    # Check if molecule matches the 2-oxo monocarboxylic acid anion pattern
    if not mol.HasSubstructMatch(oxo_acid_pattern):
         return False, "Molecule does not contain the required 2-oxo monocarboxylic acid anion substructure"
    
    # Get the matching atoms for the oxo_acid_pattern
    match_atoms = mol.GetSubstructMatches(oxo_acid_pattern)
    
    # Get the matching atoms for the carboxylic_acid_pattern
    carboxyl_atoms = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if not match_atoms:
        return False, "No match found for the 2-oxo pattern"

    # Iterate over each matching oxo-acid pattern to ensure the oxo and carboxyl are at the correct position
    for match in match_atoms:
        # Match contains the atoms in the order defined by the SMARTS pattern [CX3](=O)[CX4][CX4](=O)[OX1H0-]
        oxo_carbon = match[0]
        
        # The next carbon should have the carboxylic acid
        carboxylic_carbon_found = False
        for carbo_match in carboxyl_atoms:
          if match[2] == carbo_match[0]:
            carboxylic_carbon_found = True
            break
        if not carboxylic_carbon_found:
            return False, "Molecule does not have the carboxylate in the correct position"

    return True, "Molecule is a 2-oxo monocarboxylic acid anion"