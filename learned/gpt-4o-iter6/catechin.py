"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are characterized by a flavan-3-ol backbone, with a chroman core and hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Correct SMARTS pattern for flavan-3-ol backbone: Chroman core and possible substitutions
    flavan_3_ol_pattern = Chem.MolFromSmarts("c1cc2c(O)cc(O)c2Oc1[C@H]3COC([C@@H]3)O")
    
    if flavan_3_ol_pattern is None:
         return (None, "SMARTS pattern could not be interpreted")
    
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol (catechin) backbone found"

    # Check for hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if len(hydroxyl_matches) < 3:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxyl_matches)}"
    
    # Check stereochemistry: Look for common chiral centers assuming 2 chiral centers in flavan-3-ol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if len(chiral_centers) < 2:
        return False, "Not enough chiral centers found"

    return True, "Contains flavan-3-ol (catechin) backbone with sufficient hydroxyl groups"