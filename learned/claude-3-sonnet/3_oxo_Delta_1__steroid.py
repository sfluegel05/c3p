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

    # More flexible steroid core pattern that allows for different bond types
    # This pattern represents the four rings (A,B,C,D) of the steroid nucleus
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    )
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 3-oxo group and Delta(1) double bond in ring A
    # This pattern looks for:
    # - A double bond between C1-C2 (first two carbons)
    # - A ketone at C3 (third carbon)
    # - All in the context of the first ring
    ring_a_pattern = Chem.MolFromSmarts(
        "[#6]~1~[#6]=,:[#6]~[#6](=O)~[#6]~[#6]~[#6]~2"
    )
    
    if not mol.HasSubstructMatch(ring_a_pattern):
        return False, "Missing required 3-oxo and Delta(1) pattern"

    # Additional validation checks
    
    # Check for reasonable molecular size
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20 or num_atoms > 100:  # Steroids typically between 20-100 atoms
        return False, "Molecule size outside typical steroid range"

    # Count carbons to ensure reasonable number for steroid
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19:  # Most steroids have at least 19 carbons
        return False, "Too few carbons for a steroid"

    # Count rings to ensure reasonable number
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 4:  # Steroids must have at least 4 rings
        return False, "Too few rings for a steroid"

    return True, "Contains steroid core with 3-oxo group and Delta(1) double bond"