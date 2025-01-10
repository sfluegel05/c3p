"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDecomposition

def is_corrinoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple[bool, str]: (is_corrinoid, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cobalt atom (common in corrinoids but not required)
    has_cobalt = any(atom.GetSymbol() == 'Co' for atom in mol.GetAtoms())
    
    # Look for the basic corrin nucleus pattern:
    # Four nitrogen atoms in a macrocyclic arrangement
    n_pattern = Chem.MolFromSmarts('[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]')
    if not mol.HasSubstructMatch(n_pattern):
        return False, "Missing required nitrogen macrocycle pattern"
    
    # Count nitrogens (corrinoids typically have 4 core nitrogens)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Too few nitrogen atoms for corrinoid structure"
    
    # Look for characteristic corrin nucleus with direct C-C bond
    corrin_core = Chem.MolFromSmarts('[#7]~[#6]-[#6]~[#7]')
    if not mol.HasSubstructMatch(corrin_core):
        return False, "Missing characteristic direct C-C bond in corrin nucleus"
    
    # Check for pyrrole-like rings
    pyrrole_pattern = Chem.MolFromSmarts('[#7]1~[#6]~[#6]~[#6]~[#6]1')
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, "Missing required pyrrole-like rings"
    
    # Check for conjugated system with =C- groups
    conjugated_pattern = Chem.MolFromSmarts('[#6]=[#6]-[#6]=[#6]')
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated system with =C- groups"
    
    # Additional checks for typical corrinoid features
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() >= 4:
        return False, "Insufficient ring count for corrinoid structure"
    
    # Calculate basic molecular properties
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 300:  # Corrinoids are typically large molecules
        return False, "Molecular weight too low for corrinoid"
    
    # Success message with details
    reason = "Contains corrin nucleus with "
    reason += "4 pyrrole-like rings, conjugated system, "
    if has_cobalt:
        reason += "cobalt metal center, "
    reason += "and characteristic macrocyclic structure"
    
    return True, reason