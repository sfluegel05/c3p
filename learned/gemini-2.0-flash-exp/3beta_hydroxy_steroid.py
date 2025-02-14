"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group in the beta position on the 3rd carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid core (4 fused rings)
    # This pattern looks for the typical steroid ring system
    steroid_core_pattern = Chem.MolFromSmarts("C12[C]3[C]4[C]([C]1)[C][C]2[C]3[C]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have a steroid core structure"

    # Define the SMARTS pattern for the beta-hydroxy group at position 3
    # The [C@H] specifies that the hydrogen is below the plane (beta position)
    beta_hydroxy_pattern = Chem.MolFromSmarts("[C]12[C]3[C]4[C]([C@H](O)[C]1)[C][C]2[C]3[C]4")


    if beta_hydroxy_pattern is None:
        return False, "Invalid beta-hydroxy SMARTS pattern"


    # Check if the molecule matches the combined pattern
    matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
    
    if matches:
        return True, "Molecule matches the criteria for a 3beta-hydroxy steroid"
    else:
        return False, "Molecule does not match the steroid core with beta-hydroxyl at position 3"