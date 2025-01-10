"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
A carbohydrate derivative where one or more oxygens/hydroxy groups are replaced by sulfur or -SR
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Must contain sulfur
    if not any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "No sulfur atoms present"

    # Define sugar ring patterns more precisely
    sugar_patterns = [
        # 6-membered pyranose with proper carbon connectivity
        "[C]1[C][C][C][C]([C,O,S])1",
        # 5-membered furanose with proper carbon connectivity
        "[C]1[C][C][C]([C,O,S])1",
        # Common thiosugar patterns
        "[C]1[C][C][C][C](S[C,H])1",
        "[C]1[C][C][C](S[C,H])[C]1",
        "[C]1S[C][C][C][C]1"
    ]
    
    found_sugar_ring = False
    for pattern in sugar_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            found_sugar_ring = True
            break
            
    if not found_sugar_ring:
        return False, "No sugar ring structure found"

    # Look for direct sulfur-carbon connections characteristic of thiosugars
    thio_patterns = [
        # Direct S-C bond in ring
        "[C]1[C,O][C,O][C,O][C,O]S1",
        # Thioglycoside
        "[C]1[C,O][C,O][C,O][C,O][C]1S[C,H]",
        # Thiol group
        "[C]1[C,O][C,O][C,O][C,O][C]1S",
        # Common thiosugar linkage
        "[C]1[C,O][C,O][C,O][C,O][C]1S[C]"
    ]
    
    found_thio = False
    for pattern in thio_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            found_thio = True
            break

    # Exclude molecules where sulfur is only present as sulfate
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    if sulfate_pattern is not None:
        sulfate_matches = len(mol.GetSubstructMatches(sulfate_pattern))
        s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
        if sulfate_matches == s_count:
            return False, "Contains only sulfate groups, not thio modifications"

    # Count carbons in sugar ring
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    
    has_valid_ring = False
    for ring in ring_atoms:
        carbon_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
        if 4 <= carbon_count <= 6:  # Valid sugar ring carbon count
            has_valid_ring = True
            break
    
    if not has_valid_ring:
        return False, "No valid sugar ring carbon framework found"

    # Final classification
    if found_sugar_ring and found_thio:
        s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
        reason = "Contains sugar ring with "
        if s_count > 1:
            reason += f"{s_count} sulfur-containing modifications"
        else:
            reason += "sulfur-containing modification"
        return True, reason
    
    return False, "Structure does not match thiosugar patterns"