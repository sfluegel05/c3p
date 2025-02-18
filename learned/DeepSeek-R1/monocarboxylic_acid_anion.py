"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: CHEBI:35730 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion has a single deprotonated carboxy group (-COO-).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for exactly one carboxylate group (COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylate groups (expected 1)"
    
    # Ensure no other acidic groups (like -SO3H, -PO3H, etc.) are present
    acidic_groups = Chem.MolFromSmarts("[S,P](=O)(O)[OH]")
    if mol.HasSubstructMatch(acidic_groups):
        return False, "Contains other acidic protons"
    
    # Check if the molecule's charge is at least -1 (carboxylate contributes -1)
    charge = Chem.GetFormalCharge(mol)
    if charge > -1:
        return False, f"Formal charge {charge} is not negative enough"
    
    return True, "Single carboxylate group with no other acidic protons"