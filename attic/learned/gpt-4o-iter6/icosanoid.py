"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are derived from C20 polyunsaturated fatty acids with various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 20-carbon backbone
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbons ({c_count}), icosanoids have 20-carbon backbone"
    
    # Look for terminal carboxyl group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxyl group found"

    # Check for common functional groups
    # Hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    
    # Hydroperoxy group
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    has_hydroperoxy = mol.HasSubstructMatch(hydroperoxy_pattern)
    
    # Epoxy group
    epoxy_pattern = Chem.MolFromSmarts("[C]1OC1")
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    
    # At least one functional group should be present
    if not (has_hydroxyl or has_hydroperoxy or has_epoxy):
        return False, "Lacks functional groups typical of icosanoids"

    return True, "Contains 20-carbon backbone and typical icosanoid functional groups"