"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem contains a bicyclic ring structure with a β-lactam 
    fused to a five-membered ring and might often include a sulfur atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for carbapenem core structures
    # 4-membered ring (beta-lactam) fused to a 5-membered ring with typical substitutions
    carbapenem_patterns = [
        Chem.MolFromSmarts("C1C(=O)N2CC(C)C12"),  # Simplified for illustration
        Chem.MolFromSmarts("C1[C@H]2N(C1=O)C=C[C@H]2C")  # Bicyclic core
    ]

    matches_pattern = any(mol.HasSubstructMatch(pattern) for pattern in carbapenem_patterns)
    
    if not matches_pattern:
        return False, "Does not match any carbapenem structural core motifs"

    # Check for the presence of sulfur which is common in many carbapenems
    has_sulfur = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())

    if not has_sulfur:
        return True, "Matches carbapenem core structure but lacks typical sulfur atom"

    return True, "Matches carbapenem core structure with possible sulfur atom"

# Ensure the SMARTS patterns and additional checks comprehensively capture the desired chemical class.