"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Looks for the characteristic cyclopentanoperhydrophenanthrene skeleton
    with possible modifications.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Multiple SMARTS patterns for steroid core with variations
    steroid_patterns = [
        # Basic 6-6-6-5 fused ring system with flexibility for saturation
        "[#6]~1~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]~2~[#6]~3~4",
        # Alternative pattern allowing for some ring modifications
        "[#6]~1~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]([#6]~1)~[#6]~[#6]~[#6]~2~[#6]~3~4",
        # Pattern for possible aromatic rings
        "[#6]~1~2~[#6]~[#6]~[#6]~3~[#6]:,[#6]~[#6]:,[#6]~4~[#6]~[#6]~[#6]~[#6]([#6]~1)~[#6]~[#6]~[#6]~2~[#6]~3~4",
        # Pattern allowing for heteroatoms in rings
        "[#6,#8,#7]~1~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]([#6]~1)~[#6]~[#6]~[#6]~2~[#6]~3~4"
    ]
    
    found_skeleton = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_skeleton = True
            break
            
    if not found_skeleton:
        return False, "No steroid skeleton found"

    # Check carbon count (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Check for ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"
    
    # Check molecular weight (typical range for steroids and derivatives)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (200 <= mol_wt <= 1000):  # Expanded range to include glycosylated steroids
        return False, f"Molecular weight {mol_wt:.1f} outside typical steroid range"

    # Look for common steroid features
    features = []
    
    # Check for angular methyl groups (common but not required)
    methyl_pattern = "[CH3][C](@[#6])(@[#6])([#6])"
    if mol.HasSubstructMatch(Chem.MolFromSmarts(methyl_pattern)):
        features.append("angular methyl groups")
    
    # Common functional groups
    functional_groups = {
        "hydroxyl": "[OH]",
        "ketone": "[#6][C](=[O])[#6]",
        "ester": "[#6]C(=O)O[#6]",
        "carboxylic acid": "[#6]C(=O)[OH]",
        "double bond": "[#6]=[#6]"
    }
    
    for group_name, smarts in functional_groups.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            features.append(group_name)

    # Final classification
    reason = "Contains steroid skeleton"
    if features:
        reason += f" with {', '.join(features)}"
    
    return True, reason

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35341',
        'name': 'steroid',
        'definition': 'Any of naturally occurring compounds and synthetic analogues, '
                     'based on the cyclopenta[a]phenanthrene carbon skeleton.',
        'parents': ['CHEBI:33860']
    }
}