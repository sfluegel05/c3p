"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid consists of two or more flavonoid units connected by a single bond or atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavonoid substructure using SMARTS
    # Simplified version of a flavonoid core to match more cases
    flavonoid_pattern = Chem.MolFromSmarts(
        "[c]1[c]([c][c][c][c]1)[C](=[O])[c]2[c][c]([O])[c][c][c]2"
        )

    # Count number of flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    num_flavonoids = len(flavonoid_matches)

    if num_flavonoids < 2:
        return False, f"Found {num_flavonoids} flavonoid units, need at least 2."

    # Check if flavonoids are linked directly (or through one atom)
    # If a molecule has two or more flavonoids, and the molecule has no unconnected fragments,
    # we will assume the flavonoids are connected and therefore it's a biflavonoid.
    if len(Chem.GetMolFrags(mol)) > 1:
         return False, "Flavonoid rings are not connected, molecule is fragmented"
        
    return True, "Contains at least two flavonoid units connected by a single bond or atom"