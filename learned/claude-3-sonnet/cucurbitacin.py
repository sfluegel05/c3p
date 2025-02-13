"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:  # Lowered threshold to catch more variants
        return False, "Too few carbons for cucurbitacin"
    
    if o_count < 2:
        return False, "Too few oxygens for cucurbitacin"

    # Check for ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # More flexible tetracyclic core patterns
    core_pattern1 = Chem.MolFromSmarts("C1~C~C~2~C~C~C~3~C~4~C~C~C~3~C~2~C~1") # Basic 4-ring system
    core_pattern2 = Chem.MolFromSmarts("C1~C~C~2~C(=C)~C~C~3~C~4~C~C~C~3~C~2~C~1") # With potential double bond
    
    if not (mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2)):
        return False, "Missing characteristic ring system"

    # Check for characteristic functional groups
    
    # α,β-unsaturated ketone (common in ring A)
    unsaturated_ketone = Chem.MolFromSmarts("C=CC(=O)C")
    
    # Various oxygen-containing groups
    ketone = Chem.MolFromSmarts("C(=O)C")
    hydroxyl = Chem.MolFromSmarts("CO")
    
    # Side chain patterns (more variations)
    side_chain_patterns = [
        Chem.MolFromSmarts("CC(C)(O)CC"),
        Chem.MolFromSmarts("CC(C)(O)C=C"),
        Chem.MolFromSmarts("CC(C)(OC(=O)C)C=C"),
        Chem.MolFromSmarts("CC(O)CC=C"),
    ]
    
    has_side_chain = any(mol.HasSubstructMatch(pattern) for pattern in side_chain_patterns)
    
    # Count functional groups
    ketone_count = len(mol.GetSubstructMatches(ketone))
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl))
    
    if ketone_count + hydroxyl_count < 2:
        return False, "Insufficient oxygen-containing groups"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for cucurbitacin"

    # Count 6-membered rings
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    six_membered_rings = sum(1 for size in ring_sizes if size == 6)
    
    if six_membered_rings < 2:
        return False, "Insufficient number of 6-membered rings"

    # Scoring system
    score = 0
    score += 1 if has_side_chain else 0
    score += 1 if mol.HasSubstructMatch(unsaturated_ketone) else 0
    score += 1 if ketone_count >= 2 else 0
    score += 1 if hydroxyl_count >= 2 else 0
    score += 1 if six_membered_rings >= 3 else 0
    
    if score >= 3:
        return True, "Contains cucurbitacin ring system and characteristic functional groups"
    else:
        return False, "Missing multiple characteristic cucurbitacin features"