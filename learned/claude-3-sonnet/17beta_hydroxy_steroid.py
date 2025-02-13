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

    # More flexible steroid core pattern that allows for variations
    # This pattern matches the basic 4-ring system with more flexibility
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6,#7]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~1~[#6]~2~[#6]~3~4"
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 17β-hydroxy group
    # Matches a hydroxyl group attached to a carbon that's part of the D ring
    # The [C;R] ensures the carbon is part of a ring
    # The @@ specifies the beta stereochemistry
    oh_17_beta = Chem.MolFromSmarts('[C;R](@[#6])(@[#6])(@[#6])[OH1]')
    
    if not mol.HasSubstructMatch(oh_17_beta):
        return False, "No hydroxyl group with correct connectivity found"

    # Get all matches of the steroid core
    core_matches = mol.GetSubstructMatches(steroid_core)
    if not core_matches:
        return False, "Could not map steroid core structure"

    # Get matches for OH group
    oh_matches = mol.GetSubstructMatches(oh_17_beta)
    
    # Check for 17β-hydroxy
    found_17_beta_oh = False
    for core_match in core_matches:
        for oh_match in oh_matches:
            c_atom = mol.GetAtomWithIdx(oh_match[0])
            
            # Check if the carbon with OH is part of the D ring
            # In the steroid core pattern, the D ring carbons are typically in positions 13-17
            d_ring_positions = set(core_match[12:17])  # Indices 12-16 correspond to D ring
            if oh_match[0] in d_ring_positions:
                # Verify beta stereochemistry
                if c_atom.GetChiralTag() in [Chem.ChiralType.CHI_TETRAHEDRAL_CW, 
                                           Chem.ChiralType.CHI_TETRAHEDRAL_CCW]:
                    # Count connections to verify it's position 17
                    # Position 17 typically has 3 connections to carbons and 1 to OH
                    carbon_neighbors = sum(1 for neighbor in c_atom.GetNeighbors() 
                                        if neighbor.GetAtomicNum() == 6)
                    if carbon_neighbors == 3:
                        found_17_beta_oh = True
                        break

    if not found_17_beta_oh:
        return False, "No 17-beta hydroxyl group found"

    # Additional validation
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 18 or num_atoms > 100:  # Adjusted minimum size
        return False, "Molecule size not consistent with steroid structure"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    return True, "Contains steroid core with 17-beta hydroxyl group"