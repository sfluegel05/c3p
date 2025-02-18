"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose has six carbons and D-configuration at position 5.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for exactly 6 carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbons, found {c_count}"
    
    # Find CH2OH group (carbon with two Hs and an -OH)
    # SMARTS for -CH2-OH (not part of a ring)
    ch2oh_pattern = Chem.MolFromSmarts("[CH2][OH]")
    ch2oh_matches = mol.GetSubstructMatches(ch2oh_pattern)
    if not ch2oh_matches:
        return False, "No CH2OH group found"
    
    # Get the carbon adjacent to CH2OH (position 5)
    position5_candidates = set()
    for match in ch2oh_matches:
        ch2_atom = match[0]
        for neighbor in mol.GetAtomWithIdx(ch2_atom).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != ch2_atom:
                position5_candidates.add(neighbor.GetIdx())
    
    if not position5_candidates:
        return False, "No adjacent carbon to CH2OH found"
    
    # Check if any candidate is a chiral center with R configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers found"
    
    # Check each candidate carbon for being a chiral center with R configuration
    d_config_found = False
    for center in chiral_centers:
        atom_idx, cip = center
        if atom_idx in position5_candidates and cip == 'R':
            d_config_found = True
            break
    
    if not d_config_found:
        return False, "Position 5 does not have R configuration"
    
    # Additional check for hydroxyl groups (at least 4)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 4:
        return False, f"Expected at least 4 hydroxyls, found {hydroxyl_count}"
    
    return True, "D-configuration at position 5 (R configuration)"