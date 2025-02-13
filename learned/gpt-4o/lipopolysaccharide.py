"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Lipopolysaccharides are complex molecules consisting of a trisaccharide repeating unit
    (two heptose units and octulosonic acid) with oligosaccharide side chains and 
    3-hydroxytetradecanoic acid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool or None: True if molecule is a lipopolysaccharide, False otherwise or None if classification is not possible
        str or None: Reason for classification or None if not classified
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for key structures
    patterns = {
        "heptose": Chem.MolFromSmarts("CC(O)C(C(C(C=O)O)O)O"), # Simplified heptose pattern
        "octulosonic_acid": Chem.MolFromSmarts("O=C(O)C(C(C(C=O)O)O)CO"), # Simplified octulosonic acid pattern
        "fatty_acid": Chem.MolFromSmarts("CCCCCCCCCCCC[CH](O)C(=O)O"), # Simplified pattern for 3-hydroxytetradecanoic acid
    }

    matched_patterns = []
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.append(name)

    # Check for minimal criteria satisfaction
    if "heptose" in matched_patterns and "octulosonic_acid" in matched_patterns and "fatty_acid" in matched_patterns:
        return True, "Matches indicative motifs of a lipopolysaccharide"

    if not matched_patterns:
        return False, "No matching indicative patterns found. Complexity likely surpasses current SMILES-based method."

    return False, "Insufficient pattern match certainty for classification"