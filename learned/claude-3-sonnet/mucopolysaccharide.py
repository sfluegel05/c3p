"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating units from 
    uronic acids and glycosamines, commonly partially esterified with sulfuric acid.
    
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

    # Helper function to safely create and check SMARTS patterns
    def create_pattern(smarts):
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            print(f"Warning: Invalid SMARTS pattern: {smarts}")
            return None
        return pattern

    # Basic structural patterns
    patterns = {
        'sugar_ring': '[C;R][C;R][C;R][C;R][C;R]',  # Basic ring
        'uronic': '[C;R]C(=O)O',  # Carboxylic acid
        'amine': '[C;R][NH2,NH]',  # Amine group
        'glycosidic': '[C;R]O[C;R]',  # Glycosidic linkage
        'sulfate': 'OS(=O)(=O)O',  # Sulfate group
        'amide': '[C;R]NC(=O)',  # Amide linkage
        'hydroxy': '[C;R][OH]'  # Hydroxyl group
    }
    
    # Create patterns
    valid_patterns = {}
    for name, smarts in patterns.items():
        pattern = create_pattern(smarts)
        if pattern is not None:
            valid_patterns[name] = pattern

    # Count matches for each pattern
    matches = {}
    for name, pattern in valid_patterns.items():
        matches[name] = len(mol.GetSubstructMatches(pattern))

    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Basic element counts
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Scoring system
    score = 0
    if matches.get('sugar_ring', 0) >= 1: score += 2
    if matches.get('uronic', 0) >= 1: score += 2
    if matches.get('amine', 0) >= 1: score += 1
    if matches.get('glycosidic', 0) >= 1: score += 1
    if matches.get('sulfate', 0) >= 1: score += 1
    if matches.get('hydroxy', 0) >= 2: score += 1
    if ring_count >= 2: score += 1
    if o_count >= 6: score += 1
    if n_count >= 1: score += 1
    if mol_wt >= 500: score += 1

    # Build detailed reason
    features = []
    if matches.get('sugar_ring', 0) >= 1:
        features.append(f"{matches['sugar_ring']} sugar-like rings")
    if matches.get('uronic', 0) >= 1:
        features.append(f"{matches['uronic']} uronic acid groups")
    if matches.get('amine', 0) >= 1:
        features.append(f"{matches['amine']} amine groups")
    if matches.get('sulfate', 0) >= 1:
        features.append(f"{matches['sulfate']} sulfate groups")

    # Classification logic
    if score >= 7:
        reason = "Contains " + ", ".join(features)
        return True, reason
    
    # Rejection with specific reason
    if matches.get('sugar_ring', 0) == 0:
        return False, "No sugar-like ring structures found"
    if matches.get('uronic', 0) == 0:
        return False, "No uronic acid components found"
    if matches.get('amine', 0) == 0 and matches.get('amide', 0) == 0:
        return False, "No amine/amide components found"
    
    return False, f"Insufficient structural features (score: {score}/11)"