"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:24863 organoiodine compound
Definition: A compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound must contain at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all iodine atoms
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    
    if not iodine_atoms:
        return False, "No iodine atoms found"

    # Check if any iodine is bonded to carbon
    for iodine in iodine_atoms:
        for neighbor in iodine.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atomic number
                return True, "Contains at least one carbon-iodine bond"
    
    return False, "No carbon-iodine bonds found"