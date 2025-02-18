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
    acidic_patterns = [
        ("[CX3](=O)[OH]", "carboxylic acid (-COOH)"),
        ("[c;!$(c-[OX2])][OH]", "phenolic -OH"),  # Phenol where O is not part of a carboxylate
        ("[SX4](=O)(=O)([OH])", "sulfonic acid (-SO3H)"),
        ("[PX4](=O)([OH])([OH])", "phosphoric acid (-PO(OH)2)"),
        ("[NX3](=O)[OH]", "nitro group with -OH")
    ]
    
    for smarts, desc in acidic_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # Skip invalid patterns (unlikely if SMARTS is correct)
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains acidic group: {desc}"
    
    return True, "Single carboxylate group with no other acidic protons"