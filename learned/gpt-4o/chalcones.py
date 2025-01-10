"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone with a 1,3-diphenylpropenone structure (ArCH=CH(=O)Ar) and its derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible α,β-unsaturated ketone pattern
    alpha_beta_unsat_ketone_pattern = Chem.MolFromSmarts("C=CC(=O)[#6]")
    if not mol.HasSubstructMatch(alpha_beta_unsat_ketone_pattern):
        return False, "No α,β-unsaturated carbonyl group found"
    
    # Define aryl rings attached pattern connected to the ketone
    # Using [*] allows for variation in attachment
    aryl_ring1_pattern = Chem.MolFromSmarts("[cR1]1[cR1][cR1][cR1][cR1][cR1]1")
    aryl_ring2_pattern = Chem.MolFromSmarts("c1[cR1][cR1][cR1][cR1][cR1]1")
    
    # Check for correct connection
    connected_pattern = Chem.MolFromSmarts("[cR1]1[cR1][cR1][cR1][cR1][cR1]1C=CC(=O)[cR1]1[cR1][cR1][cR1][cR1][cR1]1")
    
    if not mol.HasSubstructMatch(connected_pattern):
        return False, "Aryl rings not properly attached to α,β-unsaturated carbonyl, or ketone functionality altered"
    
    return True, "Contains the structural motif of a chalcone with α,β-unsaturated carbonyl group and connected aryl rings"