"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the sugar ring structure with the carboxylate group
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)C(=O)[O-]")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "Glucuronic acid moiety not found"

    # Check if it is attached to an aromatic or large framework through ether/ester linkage
    attachment_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]")
    if not mol.HasSubstructMatch(attachment_pattern):
        return False, "Glucuronic acid not attached via the expected linkage"

    return True, "Contains a beta-D-glucuronic acid moiety attached to a larger structure via an ether/ester linkage"