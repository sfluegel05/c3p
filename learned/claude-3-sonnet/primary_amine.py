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
    # [NX3H2] = nitrogen with 2 hydrogens and 3 connections total
    # [#6] = any carbon atom
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][#6]")
    
    # Exclusion patterns
    alpha_amino_acid = Chem.MolFromSmarts("[NX3H2][CX4H1]([$([CX3](=[OX1])[OX2H,OX1-]),$([CX3](=[OX1])[OX2][#6])])")
    amide_pattern = Chem.MolFromSmarts("[NX3H2][CX3](=[OX1])")
    guanidine_pattern = Chem.MolFromSmarts("[NX3H2]C(=[NH2X3+,NHX2+0])[NH2X3]")
    n_oxide_pattern = Chem.MolFromSmarts("[NX3H2][#7^2]")
    hydrazine_pattern = Chem.MolFromSmarts("[NX3H2][NX3]")
    
    # Check for primary amine pattern
    if not mol.HasSubstructMatch(primary_amine_pattern):
        return False, "No primary amine groups found"
    
    # Get all matches
    amine_matches = mol.GetSubstructMatches(primary_amine_pattern)
    
    # Get exclusion matches
    alpha_amino_matches = set() if alpha_amino_acid is None else {m[0] for m in mol.GetSubstructMatches(alpha_amino_acid)}
    amide_matches = set() if amide_pattern is None else {m[0] for m in mol.GetSubstructMatches(amide_pattern)}
    guanidine_matches = set() if guanidine_pattern is None else {m[0] for m in mol.GetSubstructMatches(guanidine_pattern)}
    n_oxide_matches = set() if n_oxide_pattern is None else {m[0] for m in mol.GetSubstructMatches(n_oxide_pattern)}
    hydrazine_matches = set() if hydrazine_pattern is None else {m[0] for m in mol.GetSubstructMatches(hydrazine_pattern)}
    
    # Convert amine matches to set of nitrogen atoms
    amine_nitrogens = {match[0] for match in amine_matches}
    
    # Remove excluded nitrogens
    excluded_nitrogens = alpha_amino_matches | amide_matches | guanidine_matches | n_oxide_matches | hydrazine_matches
    valid_primary_amines = amine_nitrogens - excluded_nitrogens
    
    if not valid_primary_amines:
        return False, "No valid primary amine groups found"
        
    # Validate each remaining nitrogen
    for n_idx in valid_primary_amines:
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Must have exactly 2 hydrogens
        if n_atom.GetTotalNumHs() != 2:
            continue
            
        # Must have neutral formal charge
        if n_atom.GetFormalCharge() != 0:
            continue
            
        # Must have exactly one non-H neighbor
        heavy_neighbors = [neigh for neigh in n_atom.GetNeighbors()]
        if len(heavy_neighbors) != 1:
            continue
            
        # The neighbor must be carbon
        carbon_neighbor = heavy_neighbors[0]
        if carbon_neighbor.GetAtomicNum() != 6:
            continue
            
        # If we get here, we've found a valid primary amine
        return True, "Contains primary amine group(s)"
    
    return False, "No valid primary amine groups found"