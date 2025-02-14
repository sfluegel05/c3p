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

    # Pattern for octanoyl part of ester
    octanoyl_pattern = Chem.MolFromSmarts("CCCCCCC(=O)O[#6]") #Match oxygen with single bond
    octanoyl_matches = mol.GetSubstructMatches(octanoyl_pattern)
    
    octanoyl_pattern_ion = Chem.MolFromSmarts("CCCCCCC(=O)[O-]")
    octanoyl_matches_ion = mol.GetSubstructMatches(octanoyl_pattern_ion)

    if len(octanoyl_matches) == 0 and len(octanoyl_matches_ion) == 0:
         return False, "No octanoyl group found"

    return True, "Contains at least one octanoyl group connected to an ester"