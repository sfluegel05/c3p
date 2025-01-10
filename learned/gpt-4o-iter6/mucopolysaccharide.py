"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units from uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for uronic acid: Carboxylate connected to carbohydrate
    uronic_acid_pattern = Chem.MolFromSmarts("[O-]C(=O)C1COC(O)C(O)C1")  # Simplified sugar with carboxylate

    # SMARTS pattern for glycosamine: Sugar with amine
    glycosamine_pattern = Chem.MolFromSmarts("C(C([OH])C[OH])[NX3]")

    # Check for sulfate ester group
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O[CX4]")
    
    # Check for presence of patterns
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)

    # Check sequence alternation (Signature for mucopolysaccharides)
    if not uronic_acid_matches or not glycosamine_matches:
        return False, "No uronic acid or glycosamine units found"
    
    # Check for sufficient connectivity and possible alternation
    for ua in uronic_acid_matches:
        for ga in glycosamine_matches:
            if mol.GetBondBetweenAtoms(ua[0], ga[0]):
                if mol.HasSubstructMatch(sulfate_pattern):
                    return True, "Contains alternating uronic acids and glycosamines with sulfuric acid esterification"
                else:
                    return False, "Alternating structure found but lacks sulfate esterification"

    return False, "Failed to identify a characteristic alternating pattern of uronic acids and glycosamines"