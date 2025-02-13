"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11beta-hydroxy steroids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern that allows for variations
    # Matches the basic 6-6-6-5 ring system with any bonds
    steroid_core = Chem.MolFromSmarts("C1C(C)C2CCC3C4CCC(C4)C3C2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 11β-hydroxy group
    # Uses relative positions in the steroid skeleton and explicit stereochemistry
    # The [C@@H] ensures beta configuration of the OH group
    hydroxy_11beta = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[C@@H](O)~[#6]~[#6]~[#6]~3~[#6]~2~1")
    
    if not mol.HasSubstructMatch(hydroxy_11beta):
        return False, "No 11-beta hydroxy group found"

    # Additional validation checks
    
    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Too few rings for steroid structure"

    # Count carbons and oxygens
    num_carbons = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    num_oxygens = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    
    if num_carbons < 19 or num_carbons > 30:
        return False, f"Carbon count ({num_carbons}) outside typical steroid range (19-30)"
    
    if num_oxygens < 2:
        return False, "Too few oxygen atoms for 11-beta-hydroxy steroid"

    # Look for common steroid features
    # Most 11β-hydroxy steroids have additional oxygen-containing groups
    ketone = Chem.MolFromSmarts("C(=O)")
    has_ketone = mol.HasSubstructMatch(ketone)
    
    other_hydroxy = Chem.MolFromSmarts("CO")
    num_hydroxy = len(mol.GetSubstructMatches(other_hydroxy))
    
    if num_hydroxy < 2 and not has_ketone:
        return False, "Missing typical steroid oxygenated substituents"

    # Check for reasonable molecular weight range for steroids
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 600:
        return False, f"Molecular weight {mol_weight:.1f} outside typical range for steroids"

    return True, "Contains steroid core with 11-beta hydroxy group and appropriate substituents"