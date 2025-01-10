"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and vary 
    only in carbon side chains and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core patterns to account for different variants
    steroid_patterns = [
        # Standard steroid core
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        # Cyclopropane-containing core (for cycloartenol-type)
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1~[#6]~2",
        # Alternative core with different bond patterns
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]4~[#6]3~[#6]2~[#6]1"
    ]
    
    core_found = False
    for pattern in steroid_patterns:
        core = Chem.MolFromSmarts(pattern)
        if core is not None and mol.HasSubstructMatch(core):
            core_found = True
            break
    
    if not core_found:
        return False, "No steroid core structure found"

    # Look for 3β-hydroxyl group (common in phytosterols)
    hydroxyl_patterns = [
        "[#6][#6]1[#6][#6][#6][#6]2~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]2[#6]1[OH]",  # Basic pattern
        "[#6][C@@H]1[#6][#6][#6][#6]2~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]2[#6]1[OH]"  # Stereospecific
    ]
    
    hydroxyl_found = False
    for pattern in hydroxyl_patterns:
        hydroxyl = Chem.MolFromSmarts(pattern)
        if hydroxyl is not None and mol.HasSubstructMatch(hydroxyl):
            hydroxyl_found = True
            break
            
    if not hydroxyl_found:
        return False, "No characteristic 3-hydroxyl group found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Phytosterols typically have 27-30 carbons, but allow for variants
    if c_count < 25 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical range for phytosterols (25-35)"

    # Check for characteristic side chain patterns
    side_chain_patterns = [
        # Standard phytosterol side chain
        "[CH2,CH3][CH,C]([CH3,CH2])[CH2,CH][CH2,CH][CH3,CH2]",
        # Unsaturated side chain
        "[CH2,CH3][CH,C]([CH3,CH2])[CH,C]=[CH][CH3,CH2]",
        # Cyclopropyl-containing side chain
        "[CH2,CH3][C]([CH3,CH2])([CH2,CH])[CH2,CH][CH3,CH2]"
    ]
    
    side_chain_found = False
    for pattern in side_chain_patterns:
        side_chain = Chem.MolFromSmarts(pattern)
        if side_chain is not None and mol.HasSubstructMatch(side_chain):
            side_chain_found = True
            break
            
    if not side_chain_found:
        return False, "No characteristic phytosterol side chain found"

    # Check for characteristics that would exclude it from being a phytosterol
    allowed_elements = {1, 6, 8} # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Contains elements other than C, H, and O"

    # Calculate number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Count double bonds
    double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bonds > 5:
        return False, "Too many double bonds for a phytosterol"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 800:  # Range includes glycosylated forms
        return False, "Molecular weight outside typical range for phytosterols"

    return True, "Contains steroid core with characteristic phytosterol structure and 3-hydroxyl group"