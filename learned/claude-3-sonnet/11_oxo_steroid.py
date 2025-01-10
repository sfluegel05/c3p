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

    # Basic steroid core pattern - four fused rings (more flexible pattern)
    # This pattern looks for the basic connectivity of four fused rings
    # without being too specific about bond types
    steroid_core = Chem.MolFromSmarts("C1C2CCC3C4CCC(C4)C3C2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for ring system with 11-oxo group
    # This pattern specifically looks for the C and D rings with the 11-position ketone
    # The =O at position 11 is explicitly defined
    ring_cd_11_oxo = Chem.MolFromSmarts("[#6]1~[#6]~[#6]2~[#6](=O)~[#6]~[#6]~[#6]2~[#6]~[#6]1")
    
    if not mol.HasSubstructMatch(ring_cd_11_oxo):
        return False, "No ketone group at position 11"

    # Additional validation checks
    
    # Count carbons (steroids typically have 19-30 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (19 <= carbon_count <= 40):  # Expanded range to include larger derivatives
        return False, f"Carbon count ({carbon_count}) outside typical steroid range"

    # Count rings (steroids have at least 4 rings)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, f"Insufficient number of rings ({ring_count}) for steroid structure"

    # Count ketone groups (should have at least one)
    ketone_pattern = Chem.MolFromSmarts("[#6]=O")
    ketone_count = len(mol.GetSubstructMatches(ketone_pattern))
    if ketone_count < 1:
        return False, "No ketone groups found"

    # Check for characteristic angular methyl groups often present in steroids
    angular_methyl = Chem.MolFromSmarts("[CH3]C1([#6])[#6]~[#6]~[#6]2")
    if not mol.HasSubstructMatch(angular_methyl):
        # This is a warning but not a definitive rejection
        # Some steroids might not have the typical angular methyl groups
        pass

    # Verify the molecule has the right molecular weight range for steroids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mol_wt <= 800):  # Expanded range to include larger derivatives
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical steroid range"

    # If all checks pass, this is likely an 11-oxo steroid
    return True, "Molecule contains steroid core with ketone group at position 11"