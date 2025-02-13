"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:32952 primary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule contains a primary amine group based on its SMILES string.
    A primary amine has an NH2 group connected to a carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Primary amine pattern: NH2 connected to carbon
    # [NX3H2] = nitrogen with 2 hydrogens and 3 connections total (including H)
    # [CX4] = any carbon (can be aromatic or aliphatic)
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][#6]")
    
    # Patterns for groups that would disqualify a primary amine
    amide_pattern = Chem.MolFromSmarts("[NX3H2][CX3](=[OX1])")
    azide_pattern = Chem.MolFromSmarts("[#7]=[#7]=[#7]")
    diazo_pattern = Chem.MolFromSmarts("[#7]=[#7]")
    
    # Check if molecule has primary amine pattern
    if not mol.HasSubstructMatch(primary_amine_pattern):
        return False, "No primary amine groups found"
        
    # Get all matches
    amine_matches = mol.GetSubstructMatches(primary_amine_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern) if amide_pattern else []
    
    # Convert matches to sets of nitrogen atoms
    amine_nitrogens = {match[0] for match in amine_matches}
    amide_nitrogens = {match[0] for match in amide_matches}
    
    # Remove amide nitrogens
    true_primary_amines = amine_nitrogens - amide_nitrogens
    
    if not true_primary_amines:
        return False, "No valid primary amine groups found"
    
    # Check each nitrogen
    for n_idx in true_primary_amines:
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Must have exactly one non-H neighbor
        heavy_neighbors = [neighbor for neighbor in n_atom.GetNeighbors()]
        if len(heavy_neighbors) != 1:
            continue
            
        # The neighbor must be carbon
        if heavy_neighbors[0].GetAtomicNum() != 6:
            continue
            
        # Check formal charge
        if n_atom.GetFormalCharge() != 0:
            continue
            
        # If we get here, we've found a valid primary amine
        return True, "Contains primary amine group(s)"
            
    return False, "No valid primary amine groups found"