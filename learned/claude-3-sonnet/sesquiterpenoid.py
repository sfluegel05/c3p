"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:35189 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are terpenoids derived from sesquiterpenes (C15 skeleton),
    possibly modified by rearrangement or removal of some carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Sesquiterpenoids should have around 15 carbons (allowing some variation due to modifications)
    if c_count < 12 or c_count > 21:
        return False, f"Carbon count ({c_count}) outside typical range for sesquiterpenoids (12-21)"

    # Calculate rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    if ring_count == 0:
        return False, "No rings found - sesquiterpenoids are typically cyclic"
    
    # Most sesquiterpenoids have 1-4 rings
    if ring_count > 4:
        return False, f"Too many rings ({ring_count}) for a typical sesquiterpenoid"

    # Check molecular weight - should typically be between 200-400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for sesquiterpenoids"

    # Look for typical functional groups found in sesquiterpenoids
    functional_groups = {
        'alcohol': '[OH]',
        'ketone': '[CX3](=[OX1])[CX4]',
        'carboxylic_acid': '[CX3](=O)[OX2H1]',
        'ester': '[#6][CX3](=O)[OX2H0][#6]',
        'double_bond': '[CX3]=[CX3]'
    }
    
    found_groups = []
    for group_name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            found_groups.append(group_name)
    
    if not found_groups:
        return False, "No typical sesquiterpenoid functional groups found"

    # Look for branching methyl groups - common in terpenoids
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    
    if methyl_count == 0:
        return False, "No methyl groups found - unusual for sesquiterpenoids"

    # Calculate degree of unsaturation
    double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    if double_bonds == 0 and ring_count < 2:
        return False, "Insufficient unsaturation for sesquiterpenoid"

    # Additional check for typical carbon framework
    if c_count >= 12 and c_count <= 21 and ring_count >= 1 and ring_count <= 4:
        features = []
        features.append(f"{c_count} carbons")
        features.append(f"{ring_count} rings")
        if o_count > 0:
            features.append(f"{o_count} oxygens")
        features.append(f"{len(found_groups)} functional groups ({', '.join(found_groups)})")
        features.append(f"{methyl_count} methyl groups")
        
        return True, f"Matches sesquiterpenoid structure with {', '.join(features)}"

    return False, "Structure does not match typical sesquiterpenoid patterns"