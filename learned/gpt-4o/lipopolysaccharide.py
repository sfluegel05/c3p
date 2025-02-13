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
        "heptose": Chem.MolFromSmarts("[C@@H]([C@H](O)[C@H](CO)O)[C@@H]([C@H](O)[C=O])O"),  # Possible refined heptose pattern
        "octulosonic_acid": Chem.MolFromSmarts("O=C(O)C([C@H](C(=O)O)O)CO"),  # Refined octulosonic acid pattern
        "fatty_acid": Chem.MolFromSmarts("CCCCCCCCCCCCCCOCCC(=O)O"),  # Extended and refined 3-hydroxytetradecanoic acid pattern
    }
    
    # Initialize pattern check counts
    heptose_count = 0
    octulosonic_acid_count = 0
    fatty_acid_count = 0
    
    # Check for each pattern
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            if name == "heptose":
                heptose_count += 1
            elif name == "octulosonic_acid":
                octulosonic_acid_count += 1
            elif name == "fatty_acid":
                fatty_acid_count += 1
    
    # Verify presence and structure
    if heptose_count >= 2 and octulosonic_acid_count >= 1 and fatty_acid_count >= 1:
        return True, "Matches key motifs of a lipopolysaccharide (heptose, octulosonic acid, and fatty acid components detected)"
    
    # Determine issues with matches
    if heptose_count < 2:
        return False, f"Insufficient heptose units, found {heptose_count}"
    if octulosonic_acid_count < 1:
        return False, "Missing octulosonic acid component"
    if fatty_acid_count < 1:
        return False, "Missing fatty acid component"
    
    return False, "Complexity likely surpasses current SMILES-based method"