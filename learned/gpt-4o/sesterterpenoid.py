"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid typically has a C25 skeleton derived from a sesterterpene.
    
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
    
    # Count carbon atoms and check for a typical terpenoid structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the carbon count is around 25, which is typical for sesterterpenoids
    if c_count < 20 or c_count > 30:
        return False, "Carbon count not typical for sesterterpenoids"
    
    # Identify if the structure includes terpenoid-like features
    # Common pattern for isoprene units in terpenoids
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")  # Isoprene unit
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene substructure match found"

    # Additional checks for typical terpenoid modifications might be placed here
    # Since sesterterpenoids can have various rearrangements and functional groups

    # If the molecule passes these checks, consider it to be a potential sesterterpenoid
    return True, "Molecule matches basic terpenoid structure with expected carbon count"