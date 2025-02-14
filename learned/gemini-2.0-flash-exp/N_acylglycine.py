"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is a glycine where the nitrogen atom is acylated with another group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acylglycine, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core N-acylglycine pattern
    acyl_glycine_pattern = Chem.MolFromSmarts("NC(=O)C[CX3](=O)[OX1]")

    if not mol.HasSubstructMatch(acyl_glycine_pattern):
            return False, "Molecule does not contain the N-acylglycine core structure"

    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    
    if len(glycine_matches) < 1:
        return False, f"Molecule must have at least one glycine moiety, but found {len(glycine_matches)}"
    
    # Check that there is exactly one acylated glycine:
    acylated_matches = mol.GetSubstructMatches(acyl_glycine_pattern)
    if len(acylated_matches) !=1:
        return False, f"Molecule must have exactly one N-acylated glycine moiety, but found {len(acylated_matches)}"

    return True, "Molecule contains a glycine with the nitrogen acylated (N-acylglycine)"