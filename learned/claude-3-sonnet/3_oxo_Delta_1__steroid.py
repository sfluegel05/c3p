"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    These are steroids with a ketone at position 3 and a double bond between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings (A,B,C,D)
    # More flexible pattern that allows for variations in bond types and substitutions
    steroid_core = Chem.MolFromSmarts(
        "[#6]1~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~[#6]~[#6]~1"
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern specifically for ring A with:
    # - Double bond between positions 1-2 (Delta-1)
    # - Ketone group at position 3
    # Using explicit connectivity to ensure correct positions
    ring_a_pattern = Chem.MolFromSmarts(
        "[#6]1=[#6]-[#6](=O)-[#6]2-[#6]-[#6]-1" # Ring A with Delta-1 and 3-oxo
    )
    
    if not mol.HasSubstructMatch(ring_a_pattern):
        return False, "Missing required Delta-1 double bond and 3-oxo group in correct positions"

    # Additional validation to ensure reasonable steroid-like structure
    
    # Count rings
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 4:
        return False, "Too few rings for a steroid structure"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19:  # Most steroids have at least 19 carbons
        return False, "Too few carbons for a steroid structure"
        
    # Count oxygens (should have at least one for the 3-oxo group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "Missing oxygen atoms (required for 3-oxo group)"

    # Check for reasonable molecular weight
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 1000:  # Typical steroid weight range
        return False, "Molecular weight outside typical steroid range"

    # Verify ketone is conjugated with double bond
    conjugated_pattern = Chem.MolFromSmarts("[#6]1=[#6]-[#6](=O)")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "3-oxo group not properly conjugated with Delta-1 double bond"

    return True, "Contains steroid core with Delta-1 double bond and 3-oxo group in correct positions"