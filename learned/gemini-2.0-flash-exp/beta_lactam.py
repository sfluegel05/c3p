"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35916 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is a four-membered ring containing a nitrogen and a carbonyl, where the nitrogen is part of an amide bond within the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-lactam, focusing on the amide bond within a 4 membered ring. 
    # Also, define patterns where the N has a single bond to carbon and a double bond to a carbon.
    # This allows for the cases where the N is in an exocyclic ring.
    beta_lactam_pattern1 = Chem.MolFromSmarts("[N]1-[C](=O)-[C]-[*]1")
    beta_lactam_pattern2 = Chem.MolFromSmarts("[N]-1[C](=O)-[*]~[*]1")
    beta_lactam_pattern3 = Chem.MolFromSmarts("[N]-1~[C](=O)-[*]~[*]1")

    # Check if the molecule has the substructure
    if mol.HasSubstructMatch(beta_lactam_pattern1) or mol.HasSubstructMatch(beta_lactam_pattern2) or mol.HasSubstructMatch(beta_lactam_pattern3):
        return True, "Contains a four-membered ring with an amide bond"
    else:
        return False, "Does not contain a four-membered ring with an amide bond"