"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose has D-configuration at position 5 (R configuration in open-chain form).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all possible CH2OH groups (carbon with at least two Hs and an -OH or substituent)
    # Modified to include cases where OH is part of a substituent (e.g., phosphate)
    ch2oh_pattern = Chem.MolFromSmarts("[CH2][OX2]")
    ch2oh_matches = mol.GetSubstructMatches(ch2oh_pattern)
    if not ch2oh_matches:
        # Check for open-chain aldehyde form (terminal -CHO and CH2OH)
        aldehyde_pattern = Chem.MolFromSmarts("[CH](=O)")
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
        if not aldehyde_matches:
            return False, "No CH2OH or aldehyde group found"
        else:
            # In aldehydo form, CH2OH is at the other end
            # Find the aldehyde carbon and trace the chain
            # This part is complex; for simplicity, assume presence of CH2OH is not required
            # but focus on chiral centers
            pass
    
    # Get potential position 5 carbons (adjacent to CH2OH group)
    position5_candidates = set()
    for match in ch2oh_matches:
        ch2_atom = match[0]
        for neighbor in mol.GetAtomWithIdx(ch2_atom).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                position5_candidates.add(neighbor.GetIdx())
    
    # Also consider potential position 5 from aldehyde form
    # (This part is not fully implemented due to complexity)
    
    if not position5_candidates:
        return False, "No adjacent carbon to CH2OH found"
    
    # Check chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers found"
    
    # Check if any candidate is a chiral center with R configuration
    d_config_found = False
    for center in chiral_centers:
        atom_idx, cip = center
        if atom_idx in position5_candidates and cip == 'R':
            d_config_found = True
            break
    
    if not d_config_found:
        return False, "Position 5 does not have R configuration"
    
    return True, "D-configuration at position 5 (R configuration)"