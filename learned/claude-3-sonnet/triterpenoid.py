"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: triterpenoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Fragments

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are derived from six isoprene units and typically have 30 carbons
    in their core structure, though derivatives may have more or fewer due to modifications.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count basic ring systems
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Count sp3 carbons
    sp3_c_pattern = Chem.MolFromSmarts('[CX4]')
    sp3_c_count = len(mol.GetSubstructMatches(sp3_c_pattern))
    
    # Count methyl groups
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Look for characteristic triterpenoid core patterns
    # Steroid-like ring system (four fused rings)
    steroid_core = Chem.MolFromSmarts('C1CC2CCC3C4CCCC4CCC3C2C1')
    # Pentacyclic ring system common in many triterpenoids
    pentacyclic_core = Chem.MolFromSmarts('C1CC2CCC3C4CCCC4CCC3C2CC1')
    
    has_core = mol.HasSubstructMatch(steroid_core) or mol.HasSubstructMatch(pentacyclic_core)
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Define base criteria for triterpenoid classification
    base_criteria = (
        ring_count >= 3 and  # Allow for some ring opening/modification
        methyl_count >= 4 and  # Most triterpenoids have multiple methyl groups
        sp3_c_count >= 15 and  # Should have significant sp3 character
        mol_wt >= 400  # Minimum weight for basic triterpenoid structure
    )
    
    # Look for glycosylation pattern (for glycosidic triterpenoids)
    glycoside_pattern = Chem.MolFromSmarts('OC1OC(CO)C(O)C(O)C1')
    is_glycosylated = mol.HasSubstructMatch(glycoside_pattern)
    
    # Adjust carbon count requirements based on glycosylation
    if is_glycosylated:
        carbon_range = (35, 70)  # Broader range for glycosylated derivatives
    else:
        carbon_range = (25, 40)  # Range for basic triterpenoids and simple derivatives
        
    carbon_ok = carbon_range[0] <= c_count <= carbon_range[1]
    
    if base_criteria and (carbon_ok or is_glycosylated):
        features = []
        if has_core:
            features.append("characteristic triterpenoid core")
        if ring_count >= 4:
            features.append(f"{ring_count} rings")
        if methyl_count >= 4:
            features.append(f"{methyl_count} methyl groups")
        features.append(f"{c_count} carbons")
        features.append(f"{sp3_c_count} sp3 carbons")
        if is_glycosylated:
            features.append("glycosylated")
            
        reason = "Has " + ", ".join(features)
        return True, reason
        
    # If molecule doesn't meet criteria, explain why
    if ring_count < 3:
        return False, f"Too few rings ({ring_count}) for a triterpenoid"
    if methyl_count < 4:
        return False, f"Too few methyl groups ({methyl_count}) for a triterpenoid"
    if not carbon_ok:
        return False, f"Carbon count {c_count} outside acceptable range for triterpenoid type"
    if sp3_c_count < 15:
        return False, f"Too few sp3 carbons ({sp3_c_count}) for triterpenoid skeleton"
    if mol_wt < 400:
        return False, "Molecular weight too low for triterpenoid"
    
    return False, "Does not match typical triterpenoid characteristics"