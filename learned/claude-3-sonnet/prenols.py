"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols (CHEBI:26244)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with one or more isoprene units and a terminal hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least one hydroxyl group
    oh_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxyl group found"
    
    # Look for terminal primary alcohol (CH2-OH)
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(terminal_oh_pattern):
        return False, "No terminal primary alcohol found"
    
    # Look for isoprene unit pattern (C=C-C(C)-C)
    # This pattern can match both E and Z configurations
    isoprene_pattern = Chem.MolFromSmarts("[CH2,CH3]-[C]=[C]-[C](-[CH3])")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    if not isoprene_matches:
        return False, "No isoprene units found"
    
    # Count carbons and check if multiple of 5 (±1 for variations)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbons for a prenol"
    
    # Count double bonds (should have at least one for isoprene unit)
    double_bonds = rdMolDescriptors.CalcNumAliphaticDoubleBonds(mol)
    if double_bonds == 0:
        return False, "No double bonds found"
    
    # Check for branching pattern characteristic of isoprene units
    methyl_branches = Chem.MolFromSmarts("[CH3]-[C]=[C]")
    if not mol.HasSubstructMatch(methyl_branches):
        return False, "Missing methyl branches characteristic of isoprene units"
    
    # Additional check for phosphate groups (some prenols can be phosphorylated)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O-])")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    
    # Calculate basic properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Construct detailed reason
    reason_parts = []
    reason_parts.append(f"Contains {len(isoprene_matches)} isoprene unit(s)")
    reason_parts.append("has terminal hydroxyl group")
    if has_phosphate:
        reason_parts.append("contains phosphate group")
    reason_parts.append(f"MW: {mol_wt:.1f}")
    
    return True, "; ".join(reason_parts)