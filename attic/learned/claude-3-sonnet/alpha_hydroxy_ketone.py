"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: alpha-hydroxy ketone
A ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule contains an alpha-hydroxy ketone group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate 2D coordinates for the molecule
    AllChem.Compute2DCoords(mol)

    # SMARTS pattern for alpha-hydroxy ketone:
    # [OH]-[CH,CH2,CH3,C]-C(=O)-[#6]
    # This matches:
    # - A hydroxy group (-OH)
    # - Connected to a carbon (can be CH, CH2, CH3 or quaternary C)
    # - Connected to a ketone carbonyl (C=O)
    # - The ketone must be connected to another carbon
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts('[OH][CH,CH2,CH3,C][CX3](=[OX1])[#6]')
    
    # Alternative pattern that also catches cases where the alpha carbon is part of a ring
    alt_pattern = Chem.MolFromSmarts('[OH][C;!$(C=O)][CX3](=[OX1])[#6]')
    
    matches = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)
    alt_matches = mol.GetSubstructMatches(alt_pattern)
    
    all_matches = set(matches).union(set(alt_matches))
    
    if len(all_matches) > 0:
        return True, f"Found {len(all_matches)} alpha-hydroxy ketone group(s)"
        
    # Check if we have a ketone at all
    ketone_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[#6]')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"
        
    # Check if we have any hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy groups found"
        
    return False, "Has ketone and hydroxy groups but not in alpha position"