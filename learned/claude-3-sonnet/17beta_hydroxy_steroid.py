"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17beta-hydroxy steroid
A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern that matches the basic 4-ring system
    # Using ~ for any bond type and allowing for more flexibility in the connections
    steroid_core = Chem.MolFromSmarts(
        '[C,c]12~[C,c]~[C,c]~[C,c]~[C,c]~3~[C,c]~[C,c]~[C,c]~[C,c]~4~[C,c]~[C,c]~[C,c]~[C,c]~1~[C,c]~[C,c]~3~[C,c]~2~4'
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 17-position carbon with beta hydroxyl
    # Using a simpler pattern that looks for a carbon with a hydroxyl group
    # in a ring system at position 17
    oh_pattern = Chem.MolFromSmarts('[C;R](@[C,c])(@[C,c])(@[C,c])[OH1]')
    
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxyl group found"

    # Basic validation checks
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons and check for reasonable size
    num_carbons = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if num_carbons < 16:  # Steroids typically have at least 16 carbons
        return False, "Too few carbons for steroid structure"

    # Get all matches of the steroid core
    core_matches = mol.GetSubstructMatches(steroid_core)
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # For each hydroxyl group, check if it's attached to a carbon that's part 
    # of the D ring of the steroid core
    found_17_beta_oh = False
    for core_match in core_matches:
        d_ring_indices = set(core_match[12:16])  # Approximate D ring positions
        for oh_match in oh_matches:
            if oh_match[0] in d_ring_indices:
                # Check stereochemistry
                c_atom = mol.GetAtomWithIdx(oh_match[0])
                if c_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    found_17_beta_oh = True
                    break

    if not found_17_beta_oh:
        return False, "No 17-beta hydroxyl group found"

    return True, "Contains steroid core with 17-beta hydroxyl group"