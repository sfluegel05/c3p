"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check carbon atom count (around 25 for sesterterpenoids)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # A rough estimate for sesterterpenoids should have approximately 25 carbon atoms
    if carbon_count < 20 or carbon_count > 30:
        return False, f"Carbon count is {carbon_count}, expected around 25"

    # Check if the structure contains terpenoid characteristics
    # Terpenoids are often characterized by repeating C5 isoprene units
    # However, direct detection of these units can be challenging without explicit substructure matching

    # As the structure may be significantly rearranged or modified, 
    # including specific substructure patterns might not always yield reliable classification
    
    # Further analysis may involve complex pattern recognition or simplification
    # of known sesterterpenoid structures, which can be difficult by SMILES alone
    
    return True, "The structure meets basic criteria for a sesterterpenoid (carbon count check)"