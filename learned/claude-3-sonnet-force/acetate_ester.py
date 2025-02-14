"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:35523 acetate ester
An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for acetate group
    acetate_pattern = Chem.MolFromSmarts("CC(=O)O[*]")
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)

    # Check for ester bond
    ester_pattern = Chem.MolFromSmarts("C(=O)O[*]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Ensure the acetate group is part of the ester bond
    for acetate_match in acetate_matches:
        for ester_match in ester_matches:
            if any(atom_idx in acetate_match for atom_idx in ester_match):
                # The acetate group and ester bond overlap, indicating an acetate ester
                return True, "Contains an acetate group as part of an ester bond"

    # If no acetate ester was found
    return False, "Does not contain an acetate group as part of an ester bond"