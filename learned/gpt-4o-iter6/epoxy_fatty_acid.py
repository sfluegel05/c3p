"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide group pattern (three-membered cyclic ether: R1C1OC1R2)
    epoxide_pattern = Chem.MolFromSmarts("[#6]1-[#8]-[#6]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms to heuristically determine long chains typical of fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Fatty acids usually have at least 10 carbon atoms
        return False, f"Too few carbons for a fatty acid, found {c_count}"

    return True, "Contains epoxide group and carboxylic acid, with sufficient carbon atoms for a fatty acid"