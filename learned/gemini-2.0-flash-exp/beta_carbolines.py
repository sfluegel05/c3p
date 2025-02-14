"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is a pyrido[3,4-b]indole or its hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core beta-carboline structure using SMARTS.
    # This is a general pattern for the core, allowing for saturation and substitution.
    # The 5-membered ring has a nitrogen
    # The connection of the 5 membered nitrogen to the 6-membered ring is important
    # The connection of the indole nitrogen to the ring carbon is important
    core_pattern = Chem.MolFromSmarts('c1cc2[nH]c3c[c,n]c(c2c1)c3') 


    if mol.HasSubstructMatch(core_pattern):
      return True, "Contains the beta-carboline core structure"

    return False, "Does not contain the beta-carboline core structure"