"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:23066 carbapenem
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems have a beta-lactam ring fused to a 5-membered ring with 
    substitutions at positions 3, 4, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible core carbapenem structure patterns to match different variants
    # Pattern 1: Basic beta-lactam fused to 5-membered ring
    core_pattern1 = Chem.MolFromSmarts("[#6]1[#6][#6][#7]2[#6](=[O])[#6]12")
    
    # Pattern 2: Alternative representation with explicit double bond
    core_pattern2 = Chem.MolFromSmarts("[#6]1[#6]=[#6][#7]2[#6](=[O])[#6]12")
    
    # Pattern 3: Another common variant
    core_pattern3 = Chem.MolFromSmarts("[#6]1[#6][#6]=[#6]2[#7]1[#6]2=[O]")

    if not (mol.HasSubstructMatch(core_pattern1) or 
            mol.HasSubstructMatch(core_pattern2) or 
            mol.HasSubstructMatch(core_pattern3)):
        return False, "Missing carbapenem core structure (fused beta-lactam and 5-membered ring)"

    # Check ring sizes
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not (4 in ring_sizes and 5 in ring_sizes):
        return False, "Must contain both 4 and 5-membered rings"

    # Check for carboxylic acid/carboxylate group (common in carbapenems)
    carboxyl_patterns = [
        Chem.MolFromSmarts("C(=O)O"),  # carboxylic acid
        Chem.MolFromSmarts("C(=O)[O-]"),  # carboxylate
    ]
    has_carboxyl = any(mol.HasSubstructMatch(pat) for pat in carboxyl_patterns)
    
    # Count basic elements typically present in carbapenems
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 1:
        return False, "Must contain at least one nitrogen atom"
    
    if o_count < 2:
        return False, "Must contain at least two oxygen atoms"

    # Check molecular weight (most carbapenems are between 250-650 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for carbapenem"
    
    # Look for common substituent patterns
    substituent_patterns = [
        (Chem.MolFromSmarts("[S]"), "sulfur-containing group"),
        (Chem.MolFromSmarts("CC(O)"), "hydroxyethyl group"),
        (Chem.MolFromSmarts("C(=O)N"), "amide group")
    ]
    
    found_substituents = [desc for pat, desc in substituent_patterns if mol.HasSubstructMatch(pat)]
    
    if not has_carboxyl and not found_substituents:
        return False, "Missing characteristic substituents"

    reason = "Contains carbapenem core structure with "
    if has_carboxyl:
        reason += "carboxylic acid group "
    if found_substituents:
        reason += f"and {', '.join(found_substituents)}"
        
    return True, reason