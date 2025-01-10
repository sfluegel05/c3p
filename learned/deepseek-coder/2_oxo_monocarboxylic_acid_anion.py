"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:78289 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion is defined as an oxo monocarboxylic acid anion 
    in which the oxo group is located at the 2-position.

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

    # Check if the molecule is an anion (has a negative charge)
    negative_charges = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    if not any(charge == -1 for charge in negative_charges):
        return False, "Molecule is not an anion"

    # Check for the presence of a carboxylate group (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group found"

    # Check for the presence of an oxo group (C=O) at the 2-position relative to the carboxylate
    # The oxo group should be on the carbon adjacent to the carboxylate carbon
    oxo_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4H0][CX3](=[OX1])[O-]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        # Try a more flexible pattern that allows for different connectivity
        oxo_pattern_flexible = Chem.MolFromSmarts("[CX3](=[OX1])[CX4H0][CX3](=[OX1])")
        oxo_matches_flexible = mol.GetSubstructMatches(oxo_pattern_flexible)
        if not oxo_matches_flexible:
            return False, "No oxo group found at the 2-position relative to the carboxylate"

    # Additional check to ensure the oxo group is at the 2-position
    for match in oxo_matches:
        carboxylate_carbon = match[0]
        oxo_carbon = match[2]
        # Ensure the oxo carbon is adjacent to the carboxylate carbon
        if mol.GetBondBetweenAtoms(carboxylate_carbon, oxo_carbon) is None:
            return False, "Oxo group is not at the 2-position relative to the carboxylate"

    return True, "Contains a carboxylate group with an oxo group at the 2-position"