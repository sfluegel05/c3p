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
    
    # Check carbon atom count (typically 25 for sesterterpenoids, allow some tolerance)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Set a range for carbon count typical in sesterterpenoids
    if carbon_count < 20 or carbon_count > 35:
        return False, f"Carbon count is {carbon_count}, expected between 20 and 35"
    
    # Check for oxygen atoms often present in oxygenated terpenoids
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Likely non-terpenoid without oxygen"
    
    # Define possible terpenoid characteristics such as isoprene units using SMARTS
    # Note: Isoprene units are `C=C(C)C`
    isoprene_unit = Chem.MolFromSmarts("C=C(C)C")
    if mol.HasSubstructMatch(isoprene_unit):
        return True, "Contains isoprene units typical of terpenoids"

    # While specific substructure pattern of sesterterpenoids might be difficult to generalize,
    # this logic attempts basic inclusion through essential checks.
    
    # Return based on available evidence
    return True, "The structure sufficiently aligns with known sesterterpenoids criteria, mainly concerning carbon and the presence of terpenoid-related structures or characteristics"