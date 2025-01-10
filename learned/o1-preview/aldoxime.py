"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:16445 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as an oxime derived from an aldehyde, having the functional group R-CH=N-OH,
    where the carbon double-bonded to nitrogen is connected to only one other heavy atom (non-hydrogen),
    distinguishing it from ketoximes where the carbon is connected to two other heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the oxime functional group (C=N-OH)
    oxime_pattern = Chem.MolFromSmarts('[C]=N[OH]')
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
        return False, "No oxime functional group found"

    for match in matches:
        carbon_idx = match[0]
        nitrogen_idx = match[1]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        
        # Get neighbors of the carbon excluding the nitrogen in C=N-OH
        neighbor_atoms = [nbr for nbr in carbon.GetNeighbors() if nbr.GetIdx() != nitrogen_idx]
        
        # Count the number of heavy atom neighbors (non-hydrogen atoms)
        num_heavy_atoms = sum(1 for atom in neighbor_atoms if atom.GetAtomicNum() > 1)
        
        # Debugging statements (can be commented out in production)
        # print(f"Atom idx: {carbon_idx}, Heavy neighbors: {num_heavy_atoms}")
        
        if num_heavy_atoms == 1:
            # Carbon is attached to one other heavy atom besides nitrogen (aldoxime)
            return True, "Contains aldoxime functional group (-CH=N-OH)"
        else:
            continue  # Continue checking other matches

    return False, "No aldoxime functional group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16445',
        'name': 'aldoxime',
        'definition': 'Oximes of aldehydes RCH=NOH.',
        'parents': []
    },
}