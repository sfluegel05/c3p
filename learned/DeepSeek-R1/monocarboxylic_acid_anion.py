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
    A monocarboxylic acid anion has a single deprotonated carboxy group (-COO-) and no other acidic protons.
    
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
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups (expected 1)"
    
    # Check for other acidic groups that could donate protons (carboxylic acids, phenols, sulfonic/phosphoric acids)
    acidic_smarts = Chem.MolFromSmarts(
        "[CX3](=O)[OH] | c[OH] | [S,P](=O)(O)[OH] | [NX3](=O)[OH]"
    )
    if mol.HasSubstructMatch(acidic_smarts):
        return False, "Contains other acidic protons"
    
    # Verify formal charge is exactly -1
    charge = Chem.GetFormalCharge(mol)
    if charge != -1:
        return False, f"Formal charge {charge} does not equal -1"
    
    return True, "Single carboxylate group with no other acidic protons"