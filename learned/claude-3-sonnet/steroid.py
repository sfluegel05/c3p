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

    # Simpler SMARTS patterns for steroid core components
    core_patterns = [
        # Basic 6-6-6-5 ring connectivity with flexibility
        "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6][#6]([#6]1)[#6][#6][#6]2[#6]3",
        # Alternative pattern with more flexibility
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6][#6]([#6]1)[#6][#6][#6]2[#6]3",
        # Pattern for partially unsaturated core
        "[#6]1[#6][#6]2[#6]=,:[#6][#6]3[#6][#6][#6][#6](:[#6]1)[#6][#6][#6]2[#6]3"
    ]
    
    found_core = False
    valid_pattern = None
    for pattern in core_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            found_core = True
            valid_pattern = pattern_mol
            break
            
    if not found_core:
        return False, "No steroid core skeleton found"

    # Count rings (steroids typically have 4 or more rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Look for characteristic features
    features = []
    
    # Common substitution patterns
    patterns = {
        "angular_methyl": "[CH3][C](@[#6])(@[#6])[#6]",  # Angular methyl groups
        "hydroxyl": "[OH1]",  # Hydroxyl groups
        "ketone": "[#6]C(=O)[#6]",  # Ketone groups
        "ester": "[#6]C(=O)O[#6]",  # Ester groups
        "double_bond": "[#6]=[#6]",  # Double bonds
        "carboxyl": "[#6]C(=O)[OH1]",  # Carboxylic acid
        "alkyl_side_chain": "[CH2][CH2][CH2][CH1]([CH3])*"  # Common side chain pattern
    }

    for name, smarts in patterns.items():
        pattern_mol = Chem.MolFromSmarts(smarts)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            if name == "angular_methyl":
                features.append("angular methyl groups")
            elif name == "hydroxyl":
                features.append("hydroxyl groups")
            elif name == "ketone":
                features.append("ketone groups")
            elif name == "ester":
                features.append("ester groups")
            elif name == "double_bond":
                features.append("unsaturation")
            elif name == "carboxyl":
                features.append("carboxylic acid")
            elif name == "alkyl_side_chain":
                features.append("alkyl side chain")

    # Check number of matches to core pattern
    if valid_pattern:
        matches = len(mol.GetSubstructMatches(valid_pattern))
        if matches > 1:
            return False, "Multiple steroid-like cores found, likely not a steroid"

    # Final classification
    reason = "Contains steroid core skeleton"
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