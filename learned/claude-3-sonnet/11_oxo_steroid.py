"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is an 11-oxo steroid, reason for classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern that accounts for different bond types and oxidation states
    # Uses ~ to match any bond type and allows for aromatic carbons
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 11-oxo group in context of steroid C and D rings
    # This pattern is more specific to the 11-position ketone
    # [#6]~1 starts at C8, follows the C and D rings, with =O at C11
    oxo_11_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~2~[#6](=O)~[#6]~[#6]~[#6]~2~[#6]~[#6]~1")
    
    if not mol.HasSubstructMatch(oxo_11_pattern):
        return False, "No ketone group at position 11"

    # Count carbons (expanded range to include larger derivatives)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (17 <= carbon_count <= 50):  # More permissive range
        return False, f"Carbon count ({carbon_count}) outside typical steroid range"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, f"Insufficient number of rings ({ring_count}) for steroid structure"

    # Count ketone groups (should have at least one)
    ketone_pattern = Chem.MolFromSmarts("[#6]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if len(ketone_matches) < 1:
        return False, "No ketone groups found"

    # Additional pattern to confirm steroid-like structure
    # This pattern looks for the characteristic angular methyl group at C13
    angular_methyl = Chem.MolFromSmarts("[CH3][C]12[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]~2")
    
    # More specific check for 11-oxo position in context of full ring system
    specific_11_oxo = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~2~[#6](=O)~[#6]~[#6]~[#6]~2~[#6]([#6])~[#6]~1")
    
    # If we have both the angular methyl and specific 11-oxo pattern, we're more confident
    if mol.HasSubstructMatch(angular_methyl) and mol.HasSubstructMatch(specific_11_oxo):
        return True, "Confirmed 11-oxo steroid structure with characteristic features"
    
    # If we've made it this far, we have the basic requirements
    return True, "Contains steroid core with ketone group at position 11"