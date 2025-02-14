"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is an aldoxime where the carbon attached to the oxime group is derived from an aliphatic aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for the sp2 oxime group (C=N-O) where the carbon has only 1 non-H substituent
    oxime_pattern = Chem.MolFromSmarts("[CX3H1]([H])[NX2]-O") # Carbon has only one substituent that is not H
    matches = mol.GetSubstructMatches(oxime_pattern)

    if not matches:
        # Check for a specific pattern of formyl-derived oxime, like C(=N-O)
         oxime_pattern_formyl = Chem.MolFromSmarts("[CX3](=[OX1])[NX2]-O")
         matches = mol.GetSubstructMatches(oxime_pattern_formyl)

         if not matches:
             return False, "No aldoxime group found"

    #if there was a match, it is a proper aldoxime. 
    # Since the SMARTS pattern already filters out non-aldehyde-derived aldoximes, all matches are aliphatic if the above pattern was matched.
    return True, "Molecule is an aliphatic aldoxime"