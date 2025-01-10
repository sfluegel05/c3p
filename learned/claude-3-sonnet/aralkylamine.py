"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create SMARTS patterns
    patterns = {
        "amine": "[NX3;H2,H1,H0;!$(NC=O);!$(N=*);!$(N#*)]",  # Primary, secondary or tertiary amine
        "aromatic": "a1[a]:[a][a]:[a]1",  # Any 5 or 6-membered aromatic ring
        "aralkyl": "a[CX4][CX4,NX3;!$(NC=O)]",  # Aromatic-alkyl-N/C connection
        "amide": "[NX3][CX3](=[OX1])",  # Amide group
        "aniline": "c[NX3]"  # Direct aromatic amine
    }
    
    # Convert patterns to RDKit molecules
    compiled_patterns = {}
    for name, pattern in patterns.items():
        compiled_pattern = Chem.MolFromSmarts(pattern)
        if compiled_pattern is None:
            return False, f"Invalid SMARTS pattern for {name}"
        compiled_patterns[name] = compiled_pattern

    # Check for presence of amine
    if not mol.HasSubstructMatch(compiled_patterns["amine"]):
        return False, "No amine group found"
    
    # Check for presence of aromatic ring
    if not mol.HasSubstructMatch(compiled_patterns["aromatic"]):
        return False, "No aromatic ring found"
    
    # Get all amine nitrogens
    amine_matches = mol.GetSubstructMatches(compiled_patterns["amine"])
    if not amine_matches:
        return False, "No valid amine groups found"
    
    # Check for amides and exclude them
    amide_matches = mol.GetSubstructMatches(compiled_patterns["amide"])
    amide_nitrogens = set(match[0] for match in amide_matches)
    
    # Check for anilines and exclude them
    aniline_matches = mol.GetSubstructMatches(compiled_patterns["aniline"])
    aniline_nitrogens = set(match[1] for match in aniline_matches)
    
    # Get valid amine nitrogens (not part of amides or anilines)
    valid_amine_nitrogens = set(match[0] for match in amine_matches) - amide_nitrogens - aniline_nitrogens
    
    if not valid_amine_nitrogens:
        return False, "Contains only amides or anilines, no valid alkylamines"
    
    # Check for aralkyl pattern
    aralkyl_matches = mol.GetSubstructMatches(compiled_patterns["aralkyl"])
    if not aralkyl_matches:
        return False, "No alkyl linker between aromatic ring and amine/carbon"
        
    # Verify that at least one aralkyl pattern connects to a valid amine
    valid_connection = False
    for match in aralkyl_matches:
        if match[-1] in valid_amine_nitrogens:
            valid_connection = True
            break
            
    if not valid_connection:
        return False, "No valid connection between aromatic ring and amine via alkyl linker"
        
    return True, "Contains amine group connected to aromatic ring via alkyl linker"