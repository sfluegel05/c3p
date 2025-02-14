"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a 4'-hydroxyflavanone, False otherwise.
         str:  Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the flavanone core (benzopyranone system)
    flavanone_core_pattern = Chem.MolFromSmarts("[cH1]1[cH1][cH1][c]([O][c]2[cH1][cH1][c]([C](=O)[CH2]2)[cH1]1)")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Molecule does not contain the flavanone core."
    
    # SMARTS pattern for the 4'-hydroxy substitution, using the para position pattern
    hydroxy_4_pattern = Chem.MolFromSmarts("[cH1]1[cH1][cH1][c](O)[cH1][cH1]1")

    # Check for the 4'-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_4_pattern):
        return False, "Molecule does not have a hydroxyl substituent at the 4' position of the phenyl ring"


    return True, "Molecule contains a flavanone core with a hydroxy substituent at the 4' position of the phenyl ring"