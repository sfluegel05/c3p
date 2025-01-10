"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone with a 1,3-diphenylpropenone structure (ArCH=CH(=O)Ar).

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

    # Look for the α,β-unsaturated carbonyl group: C=C-C=O
    alpha_beta_unsat_ketone_pattern = Chem.MolFromSmarts("C=CC(=O)")
    if not mol.HasSubstructMatch(alpha_beta_unsat_ketone_pattern):
        return False, "No α,β-unsaturated carbonyl group found"
    
    # Look for two phenyl or aryl rings connected through the propenone
    aryl_pattern = Chem.MolFromSmarts("c1ccccc1")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern)
    if len(aryl_matches) < 2:
        return False, "Less than two aryl rings found"
    
    # Ensure that both aryl rings are appropriately connected to the propenone chain
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("C=CC(=O)c1ccccc1")):
        return False, "Aryl rings not properly attached to α,β-unsaturated carbonyl"

    return True, "Contains the structural motif of a chalcone with α,β-unsaturated carbonyl group and connected aryl rings"