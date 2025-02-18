"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol has a hydroxyl group attached to a saturated carbon with either three hydrogens or one other carbon and two hydrogens.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens and sanitize to ensure accurate hydrogen count
    mol = Chem.AddHs(mol)
    
    # Find all hydroxyl groups (O with exactly one H)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"
    
    for oh_match in hydroxyl_matches:
        o_idx = oh_match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Verify hydroxyl group (O with exactly one H)
        if o_atom.GetTotalNumHs() != 1:
            continue
        
        # Get attached carbon
        carbon_neighbors = [n for n in o_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # Alcohol must be attached to exactly one carbon
        c_atom = carbon_neighbors[0]
        
        # Check carbon is saturated (all single bonds)
        for bond in c_atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                break
        else:
            # Count hydrogens and adjacent carbons on the carbon (excluding hydroxyl O)
            h_count = c_atom.GetTotalNumHs()
            carbon_adjacent = 0
            other_adjacent = 0
            for neighbor in c_atom.GetNeighbors():
                if neighbor.GetIdx() == o_idx:
                    continue  # Skip hydroxyl oxygen
                if neighbor.GetAtomicNum() == 6:
                    carbon_adjacent += 1
                elif neighbor.GetAtomicNum() != 1:  # Non-carbon, non-hydrogen atoms
                    other_adjacent += 1
            
            # Disallow any non-carbon/hydrogen attachments
            if other_adjacent > 0:
                continue
            
            # Check primary alcohol conditions
            if (h_count == 3 and carbon_adjacent == 0) or (h_count == 2 and carbon_adjacent == 1):
                return True, "Hydroxyl group attached to a primary carbon"
    
    return False, "No primary alcohol group detected"