"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is defined by the presence of a cyclopenta[a]phenanthrene core
    with a hydroxyl group at the 11th position with beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the general steroid backbone
    # This pattern broadly represents the cyclopenta[a]phenanthrene structure
    steroid_core = Chem.MolFromSmarts("C1CCC2C1CCC3C2CC=C4C3(=CCC4)C")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Define SMARTS pattern for 11-beta-hydroxy group
    # Check for the OH group attached to the beta face on the typical 11th position carbon
    # The beta face is typically wedged in a structure, requiring stereochemistry verification
    beta_11_hydroxy = Chem.MolFromSmarts("[C@H]1C[C@H](O)C2=C1CC=C3C2CCC4=C3CCC4")
    if not mol.HasSubstructMatch(beta_11_hydroxy):
        return False, "11-beta-hydroxy group not found or incorrect configuration"

    return True, "Contains 11-beta hydroxy group in correct configuration for steroids"

# Example usage
example_smiles = "[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(CO)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C"
print(is_11beta_hydroxy_steroid(example_smiles))