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
    Looks for the characteristic cyclopenta[a]phenanthrene skeleton and its variations.
    
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

    # More flexible SMARTS patterns for steroid core
    core_patterns = [
        # Basic steroid core (saturated)
        "[C,c]1~[C,c]~[C,c]2~[C,c]~[C,c]~[C,c]3~[C,c]~[C,c]~[C,c]~[C,c]~1~[C,c]~[C,c]~[C,c]2~[C,c]3",
        # Variant allowing for ring modifications
        "[C,c]1~[C,c]~[C,c]2~[C,c]~[C,c]~[C,c]3~[C,c]~[C,c]~[C,c]~[C,c]~1~[C,c]~[C,c]2~[C,c]3",
        # More flexible pattern for modified steroids
        "[C,c]1~2~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~3~[C,c]~[C,c]~[C,c]~[C,c]~1~[C,c]~2~[C,c]3",
        # Pattern for heavily modified cores
        "[C,c]1~2~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~3~[C,c]~[C,c]~[C,c]~1~[C,c]~2~[C,c]3"
    ]

    # Check for steroid core
    found_core = False
    for pattern in core_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            found_core = True
            break

    if not found_core:
        return False, "No steroid-like ring system found"

    # Count rings
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, "Insufficient number of rings"

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
        "alkyl_side": "[CH2,CH1][CH2,CH1][CH2,CH1,CH3]"  # Side chain
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
            elif name == "alkyl_side":
                features.append("alkyl substitution")

    # Additional checks
    n_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_arom_rings > 1:
        return False, "Too many aromatic rings for typical steroid"

    # Calculate fraction of sp3 carbons (steroids are typically highly saturated)
    n_sp3 = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH4,CH3,CH2,CH1]")))
    if carbon_count > 0:
        sp3_fraction = n_sp3 / carbon_count
        if sp3_fraction < 0.3:  # Allow some unsaturation but not too much
            return False, "Too unsaturated for typical steroid structure"

    # Final classification
    reason = "Contains steroid ring system"
    if features:
        reason += f" with {', '.join(features)}"
    if n_rings > 4:
        reason += f" and {n_rings} rings"

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