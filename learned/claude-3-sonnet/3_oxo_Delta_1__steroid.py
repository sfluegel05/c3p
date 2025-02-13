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

    # Basic steroid core pattern - four fused rings with more flexibility
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~[#6]~[#6]~1"
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for ring A with Delta-1 double bond and 3-oxo group
    # More specific pattern to ensure correct positioning
    ring_a_pattern = Chem.MolFromSmarts(
        "[#6]1=,:[#6]-[#6](=[O])-[#6]2-[#6]-[#6]" # Ring A pattern
    )
    
    if not mol.HasSubstructMatch(ring_a_pattern):
        return False, "Missing required Delta-1 double bond and 3-oxo group pattern"

    # Additional check for specific Delta-1 double bond
    delta_1_pattern = Chem.MolFromSmarts("[#6]1=,:[#6]-[#6](=O)")
    if not mol.HasSubstructMatch(delta_1_pattern):
        return False, "Missing Delta-1 double bond"

    # Additional validation to ensure reasonable steroid-like structure
    
    # Count carbons (most steroids have 19-30 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:  # Allow some flexibility
        return False, "Too few carbons for a steroid structure"

    # Count rings (steroids typically have 4 or more rings)
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 4:
        return False, "Too few rings for a steroid structure"

    # Check for ketone oxygen
    ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    # Verify molecular weight is in reasonable range for steroids
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 200 or mol_weight > 1000:
        return False, "Molecular weight outside typical steroid range"

    # Check for conjugation between Delta-1 and 3-oxo
    conjugated_enone = Chem.MolFromSmarts("[#6]1=,:[#6]-[#6](=O)-[#6]")
    if not mol.HasSubstructMatch(conjugated_enone):
        return False, "Missing conjugated enone system"

    return True, "Contains steroid core with Delta-1 double bond conjugated to 3-oxo group"