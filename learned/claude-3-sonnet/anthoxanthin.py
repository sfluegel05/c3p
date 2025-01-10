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

    # Basic flavonoid core (more general pattern)
    # Matches the basic C6-C3-C6 skeleton of flavonoids
    flavonoid_core = Chem.MolFromSmarts("c1c(c(=O)c2c(o1)cc[c,n]c2)c3ccccc3")
    
    # Alternative flavonoid core patterns
    flavone_core = Chem.MolFromSmarts("O=C1C=C(Oc2ccccc12)c3ccccc3")
    flavonol_core = Chem.MolFromSmarts("O=C1C(O)=C(Oc2ccccc12)c3ccccc3")
    
    # More general benzopyrone core
    benzopyrone_core = Chem.MolFromSmarts("O=C1CCOc2c1cccc2")
    
    # Check for presence of any core structure
    core_patterns = [p for p in [flavonoid_core, flavone_core, flavonol_core, benzopyrone_core] if p is not None]
    has_core = any(mol.HasSubstructMatch(pattern) for pattern in core_patterns)
    
    if not has_core:
        return False, "Missing basic flavonoid/benzopyrone core structure"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:  # Minimum oxygens needed for basic structure
        return False, "Insufficient oxygen atoms for anthoxanthin"

    # Look for characteristic groups
    patterns = {
        'hydroxyl': Chem.MolFromSmarts("[OH]"),
        'methoxy': Chem.MolFromSmarts("OC"),
        'ketone': Chem.MolFromSmarts("C(=O)"),
        'ether': Chem.MolFromSmarts("COC"),
        'glycoside': Chem.MolFromSmarts("OC1OCC(O)C(O)C1")
    }
    
    matches = {}
    for name, pattern in patterns.items():
        if pattern:
            matches[name] = len(mol.GetSubstructMatches(pattern))

    # Must have ketone group (part of core structure)
    if matches.get('ketone', 0) < 1:
        return False, "Missing required ketone group"

    # Must have some oxygenated substituents
    total_oxy_groups = (matches.get('hydroxyl', 0) + 
                       matches.get('methoxy', 0) + 
                       (1 if matches.get('glycoside', 0) > 0 else 0))
    
    if total_oxy_groups < 1:
        return False, "Missing required oxygenated substituents"

    # Check for aromatic character
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 6:
        return False, "Insufficient aromatic character"

    # Count carbons (relaxed requirement)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:  # Typical flavonoid carbon count
        return False, "Insufficient carbon atoms for anthoxanthin structure"

    # Classification reason
    reason_parts = []
    if matches.get('hydroxyl', 0) > 0:
        reason_parts.append("hydroxyl groups")
    if matches.get('methoxy', 0) > 0:
        reason_parts.append("methoxy groups")
    if matches.get('glycoside', 0) > 0:
        reason_parts.append("glycoside moiety")
    
    reason = "Contains flavonoid core with " + ", ".join(reason_parts)
    
    return True, reason