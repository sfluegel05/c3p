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
        # Any aromatic ring (including heterocycles)
        "aromatic": "[a;r5,r6]1[a;r5,r6][a;r5,r6][a;r5,r6][a;r5,r6]1",
        
        # Amine connected to aliphatic carbon (not amide, not imine, not aniline)
        "alkylamine": "[NX3;H2,H1,H0;!$(NC=O);!$(N=*);!$(N#*);!$(Na);!$(N[a])]",
        
        # Aromatic-alkyl-amine connection pattern (flexible chain length)
        "aralkyl": "[a;r5,r6]~[C;!$(C=O);!$(C=N);!$(C#N)]~[*]~[NX3;!$(NC=O);!$(N=*);!$(N#*);!$(Na)]",
        
        # Patterns to exclude
        "amide": "[NX3][CX3](=[OX1])",
        "imine": "[NX2]=[CX3]",
        "aniline": "[a][NX3]",
        "guanidine": "[NX3][CX3](=[NX2])[NX3]"
    }
    
    # Convert patterns to RDKit molecules
    compiled_patterns = {}
    for name, pattern in patterns.items():
        compiled_pattern = Chem.MolFromSmarts(pattern)
        if compiled_pattern is None:
            return False, f"Invalid SMARTS pattern for {name}"
        compiled_patterns[name] = compiled_pattern

    # Check for presence of aromatic ring
    if not mol.HasSubstructMatch(compiled_patterns["aromatic"]):
        return False, "No aromatic ring found"

    # Check for presence of alkylamine
    if not mol.HasSubstructMatch(compiled_patterns["alkylamine"]):
        return False, "No alkylamine group found"

    # Get all nitrogen atoms
    nitrogens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    
    # Exclude nitrogens that are part of excluded patterns
    excluded_nitrogens = set()
    for pattern_name in ["amide", "imine", "aniline", "guanidine"]:
        matches = mol.GetSubstructMatches(compiled_patterns[pattern_name])
        for match in matches:
            excluded_nitrogens.update([idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7])
    
    valid_nitrogens = set(nitrogens) - excluded_nitrogens
    if not valid_nitrogens:
        return False, "No valid alkylamine groups found"

    # Check for aralkyl connection
    aralkyl_matches = mol.GetSubstructMatches(compiled_patterns["aralkyl"])
    if not aralkyl_matches:
        return False, "No connection between aromatic ring and amine via alkyl linker"

    # Verify that at least one aralkyl pattern connects to a valid amine
    for match in aralkyl_matches:
        amine_idx = match[-1]
        if amine_idx in valid_nitrogens:
            return True, "Contains amine group connected to aromatic ring via alkyl linker"

    return False, "No valid aralkylamine structure found"