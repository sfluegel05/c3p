"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments that include flavones and flavonols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Multiple core patterns to catch different flavonoid variations
    core_patterns = [
        # Basic flavone/flavonol core (more flexible)
        "O=C1CC(c2ccccc2)Oc2ccccc12",
        # Alternative core with different bond orders
        "O=C1C=C(c2ccccc2)Oc2ccccc12",
        # More general pattern for modified cores
        "O=C1C(=C)C(c2ccccc2)Oc2ccccc12",
        # Pattern for isoflavones
        "O=C1C(c2ccccc2)=COc2ccccc12",
        # Pattern catching benzopyrone core
        "O=C1CCOc2c1cccc2",
        # Pattern for modified flavonoids
        "O=C1C(O)C(c2ccccc2)Oc2ccccc12"
    ]

    has_core = False
    for pattern in core_patterns:
        core = Chem.MolFromSmarts(pattern)
        if core and mol.HasSubstructMatch(core):
            has_core = True
            break

    if not has_core:
        return False, "Missing flavonoid core structure"

    # Count key atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if o_count < 2:
        return False, "Insufficient oxygen atoms for anthoxanthin"
    
    if c_count < 15:
        return False, "Insufficient carbon atoms for anthoxanthin structure"

    # Look for characteristic groups with more flexible patterns
    patterns = {
        'hydroxyl': '[OH]',
        'methoxy': 'CO[c,C]',
        'ketone': 'C(=O)',
        'ether': 'COC',
        # More flexible glycoside pattern
        'glycoside': '[CH1,CH2]O[CH]1O[CH][CH][CH][CH][CH]1',
        # Pattern for prenyl/geranyl groups
        'prenyl': 'CC(C)=CC',
        # Pattern for aromatic rings
        'aromatic': 'c1ccccc1'
    }
    
    matches = {}
    for name, pattern in patterns.items():
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            matches[name] = len(mol.GetSubstructMatches(patt))

    # Must have ketone group and some oxygenated substituents
    if matches.get('ketone', 0) < 1:
        return False, "Missing required ketone group"

    # Check for characteristic substitution
    total_oxy_groups = (matches.get('hydroxyl', 0) + 
                       matches.get('methoxy', 0) + 
                       matches.get('glycoside', 0))
    
    if total_oxy_groups < 1:
        return False, "Missing required oxygenated substituents"

    # Check for aromatic character
    if matches.get('aromatic', 0) < 1:
        return False, "Missing required aromatic rings"

    # Build classification reason
    reason_parts = []
    if matches.get('hydroxyl', 0) > 0:
        reason_parts.append("hydroxyl groups")
    if matches.get('methoxy', 0) > 0:
        reason_parts.append("methoxy groups")
    if matches.get('glycoside', 0) > 0:
        reason_parts.append("glycosidic substituents")
    if matches.get('prenyl', 0) > 0:
        reason_parts.append("prenyl/geranyl groups")
    
    base_reason = "Contains flavonoid core structure"
    if reason_parts:
        base_reason += " with " + ", ".join(reason_parts)
    
    return True, base_reason