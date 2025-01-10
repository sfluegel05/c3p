"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxyl group at position 3 in the alpha orientation.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid nucleus (tetracyclic fused ring system)
    # Define SMARTS pattern for steroid nucleus
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3CCCC4')  # Simplified steroid core
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"
    
    # Define SMARTS pattern for 3alpha-hydroxy group
    # Assuming the 3-position is part of the first ring in the steroid nucleus
    hydroxyl_pattern = Chem.MolFromSmarts('[C@H](O)[C@@H]1CC[C@H]2[C@@H](C1)CC[C@]3(C)C2')  # Simplified pattern
    matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not matches:
        # Try matching with possible stereochemistry variations
        hydroxyl_pattern = Chem.MolFromSmarts('[C@H](O)[C@H]1CC[C@H]2[C@@H](C1)CC[C@]3(C)C2')
        matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if not matches:
            return False, "No 3alpha-hydroxy group found in steroid nucleus"
    
    # Check stereochemistry at C3
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    c3_atom_idx = None
    for match in matches:
        # Assuming the first atom in the match is C3
        c3_idx = match[0]
        c3_atom_idx = c3_idx
        break  # Consider the first match
    
    if c3_atom_idx is None:
        return False, "Unable to identify C3 position"
    
    # Get stereochemistry at C3
    for idx, chirality in chiral_centers:
        if idx == c3_atom_idx:
            if chirality == 'S':
                return True, "Molecule is a 3alpha-hydroxy steroid"
            else:
                return False, "Hydroxyl group at C3 is not in alpha position"
    return False, "Chiral center at C3 not found or stereochemistry not assigned"