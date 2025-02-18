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
    
    # Find all hydroxyl groups (O with H attached)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"
    
    for oh_match in hydroxyl_matches:
        o_idx = oh_match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Check if the oxygen is part of a hydroxyl group (bonded to exactly one H)
        if o_atom.GetTotalNumHs() != 1:
            continue
        
        # Get the carbon attached to the hydroxyl oxygen
        neighbors = o_atom.GetNeighbors()
        if len(neighbors) != 1:
            continue  # oxygen should be bonded to one carbon
        c_atom = neighbors[0]
        
        # Check if the carbon is saturated (all single bonds)
        for bond in c_atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                break
        else:
            # Count hydrogens and adjacent carbons on the carbon
            h_count = c_atom.GetTotalNumHs()
            carbon_neighbors = 0
            other_neighbors = 0
            for neighbor in c_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    carbon_neighbors += 1
                elif neighbor.GetAtomicNum() != 1:  # exclude hydrogens
                    other_neighbors += 1
            
            # Check for any non-carbon, non-hydrogen neighbors
            if other_neighbors > 0:
                continue
            
            # Check primary alcohol conditions
            if (h_count == 3 and carbon_neighbors == 0) or (h_count == 2 and carbon_neighbors == 1):
                return True, "Hydroxyl group attached to a primary carbon"
    
    return False, "No primary alcohol group detected"