"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the carboxylic acid component is octanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern to match octanoyl connected to ester bond
    # Match CCCCCCCC(=O)O connected to a C
    octanoyl_ester_pattern = Chem.MolFromSmarts("CCCCCCC(=O)O[#6]")
    octanoyl_ester_matches = mol.GetSubstructMatches(octanoyl_ester_pattern)

    # Match ionized version CCCCCCCC(=O)[O-] connected to a C
    octanoyl_ester_pattern_ion = Chem.MolFromSmarts("CCCCCCC(=O)[O-][#6]")
    octanoyl_ester_matches_ion = mol.GetSubstructMatches(octanoyl_ester_pattern_ion)


    if len(octanoyl_ester_matches) == 0 and len(octanoyl_ester_matches_ion) == 0:
        return False, "No octanoyl ester group found"

    return True, "Contains at least one octanoyl group connected to an ester"