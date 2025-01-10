"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:26154 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a furan ring attached to a steroid skeleton
    furan_steroid_pattern = Chem.MolFromSmarts("[C;H1,H2]1[C;H1,H2][C;H1,H2][C;H1,H2][O;H1]1.[C;H1,H2]2[C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;