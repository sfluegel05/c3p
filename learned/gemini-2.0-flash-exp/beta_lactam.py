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

    # Define the SMARTS pattern for beta-lactam
    beta_lactam_pattern = Chem.MolFromSmarts("[NX2;R1][CX3;R1](=[OX1])[R1]")

    # Check if the molecule has the substructure
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered ring with an amide bond"
    else:
        return False, "Does not contain a four-membered ring with an amide bond"