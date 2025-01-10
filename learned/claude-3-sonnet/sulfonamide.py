"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: sulfonamide
Definition: An amide of a sulfonic acid RS(=O)2NR'2
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide contains the RS(=O)2NR'2 group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonamide pattern: R-S(=O)(=O)-N
    # Note: We use recursive SMARTS to ensure S is connected to at least one carbon
    sulfonamide_pattern = Chem.MolFromSmarts('[C,c]S(=O)(=O)[N]')
    
    if not mol.HasSubstructMatch(sulfonamide_pattern):
        return False, "No sulfonamide group (RS(=O)2N) found"
    
    # Count number of sulfonamide groups
    matches = mol.GetSubstructMatches(sulfonamide_pattern)
    num_sulfonamide = len(matches)
    
    # Verify that sulfur has correct oxidation state and connectivity
    for match in matches:
        s_idx = match[1]  # Index of sulfur atom in the match
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Check sulfur valence and oxidation state
        if s_atom.GetTotalValence() != 6:
            return False, f"Sulfur atom has incorrect valence {s_atom.GetTotalValence()}, should be 6"
            
        # Check number of oxygen atoms connected to sulfur (should be 2)
        oxygen_count = len([n for n in s_atom.GetNeighbors() 
                          if n.GetAtomicNum() == 8 and n.GetTotalValence() == 2])
        if oxygen_count != 2:
            return False, f"Sulfur atom has {oxygen_count} double-bonded oxygens, should be 2"
    
    # Additional validation to exclude sulfonic acids and esters
    sulfonic_acid_pattern = Chem.MolFromSmarts('S(=O)(=O)O')
    sulfonic_ester_pattern = Chem.MolFromSmarts('S(=O)(=O)OC')
    
    if mol.HasSubstructMatch(sulfonic_acid_pattern):
        sulfonic_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
        # Check if all S(=O)2O groups are part of sulfonamides
        if len(sulfonic_matches) > num_sulfonamide:
            return False, "Contains sulfonic acid group"
            
    if mol.HasSubstructMatch(sulfonic_ester_pattern):
        ester_matches = mol.GetSubstructMatches(sulfonic_ester_pattern)
        # Check if all S(=O)2OC groups are part of sulfonamides
        if len(ester_matches) > num_sulfonamide:
            return False, "Contains sulfonic ester group"
    
    return True, f"Contains {num_sulfonamide} sulfonamide group(s) (RS(=O)2N)"