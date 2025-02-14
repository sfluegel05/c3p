"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Define a more general flavonoid substructure using SMARTS
    # Looking for the core benzopyran ring system
    flavonoid_pattern = Chem.MolFromSmarts(
        "c1cc(Oc2ccccc2C(=O))cc1"
    )


    # Count number of flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    num_flavonoids = len(flavonoid_matches)

    if num_flavonoids < 2:
        return False, f"Found {num_flavonoids} flavonoid units, need at least 2."

    # Check if molecule is a single fragment (to ensure connection), no more detailed analysis
    if len(Chem.GetMolFrags(mol)) > 1:
        return False, "Flavonoid rings are not connected, molecule is fragmented"
    

    return True, "Contains at least two flavonoid units connected by a single bond or atom"