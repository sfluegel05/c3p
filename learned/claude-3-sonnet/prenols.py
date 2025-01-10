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
    
    # Check for terminal primary alcohol (CH2-OH)
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(terminal_oh_pattern):
        return False, "No terminal primary alcohol found"
    
    # Count rings - prenols should be mostly linear
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 1:
        return False, f"Too many rings ({ring_count}) for a typical prenol"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbons for a prenol"
    
    # Look for isoprene units with various patterns
    isoprene_patterns = [
        # Regular isoprene unit with head-to-tail connection
        Chem.MolFromSmarts("CC(=C)CC"),
        # Trans isoprene unit
        Chem.MolFromSmarts("C/C=C(/C)CC"),
        # Cis isoprene unit
        Chem.MolFromSmarts("C\C=C(\C)CC"),
        # Saturated isoprene unit
        Chem.MolFromSmarts("CC(C)CCC"),
    ]
    
    total_isoprene_matches = 0
    for pattern in isoprene_patterns:
        if pattern is not None:
            matches = len(mol.GetSubstructMatches(pattern))
            total_isoprene_matches += matches
    
    if total_isoprene_matches == 0:
        return False, "No isoprene units found"
        
    # Check for proper carbon count multiple of 5 (allowing some deviation)
    if c_count % 5 > 2:
        return False, f"Carbon count {c_count} not consistent with isoprene units"
    
    # Look for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if double_bond_pattern:
        double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    else:
        double_bonds = 0
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for reasonable molecular weight range for prenols
    if mol_wt < 70:  # Smallest prenol (C5H10O) would be about 86
        return False, "Molecular weight too low for prenol"
        
    # Check for excessive heteroatoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if n_count > 1 or s_count > 0:
        return False, "Too many heteroatoms for a prenol"
        
    # Additional check for phosphate groups (some prenols can be phosphorylated)
    phosphate_patterns = [
        Chem.MolFromSmarts("[P](=[O])([O-])([O-])=O"),  # Diphosphate
        Chem.MolFromSmarts("[P](=[O])([O-])[O]"),       # Phosphate
    ]
    has_phosphate = any(pattern is not None and mol.HasSubstructMatch(pattern) 
                       for pattern in phosphate_patterns)
    
    # Construct detailed reason
    reason_parts = []
    reason_parts.append(f"Contains {total_isoprene_matches} isoprene unit(s)")
    reason_parts.append("has terminal hydroxyl group")
    reason_parts.append(f"has {double_bonds} double bond(s)")
    if has_phosphate:
        reason_parts.append("contains phosphate group")
    reason_parts.append(f"MW: {mol_wt:.1f}")
    
    return True, "; ".join(reason_parts)