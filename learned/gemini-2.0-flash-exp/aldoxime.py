"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core aldoxime substructure using SMARTS
    oxime_pattern = Chem.MolFromSmarts('[CH][C]=[N][O]') # C=N-O with C attached to at least 1 H
    aldehyde_carbon = Chem.MolFromSmarts('[CH1]')

    # Check if the core substructure and the aldehyde carbon are present
    if not mol.HasSubstructMatch(oxime_pattern):
            return False, "Molecule does not contain the C=N-O oxime group"
    
    # Find matches for the [CH][C]=[N][O] pattern, and ensure that the [CH] part is connected to one H
    matches = mol.GetSubstructMatches(oxime_pattern)
    for match in matches:
      c_atom = mol.GetAtomWithIdx(match[0])
      if c_atom.HasSubstructMatch(aldehyde_carbon):
          return True, "Molecule contains the R-CH=NOH aldoxime substructure."

    return False, "Molecule does not contain the correct R-CH=NOH aldoxime substructure."