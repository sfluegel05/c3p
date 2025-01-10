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
        bool: True if the molecule is an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the epoxide group pattern (three-membered cyclic ether: C1OC1)
    epoxide_pattern = Chem.MolFromSmarts("[C]1-[O]-[C]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms to heuristically determine a long chain typical of fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:  # Adjusted heuristic for fatty acids
        return False, f"Too few carbons for a fatty acid, found {c_count}"

    # Identify patterns that should not be present, such as linked epoxide rings or non-linear chains
    unwanted_patterns = [
        Chem.MolFromSmarts("[C]1-[O]-[C]2-[C]1-[O]-[C]2"),  # Connected epoxide rings
        Chem.MolFromSmarts("[C](=[O])O[C]")  # Ester linkages from carboxylate
    ]

    for pattern in unwanted_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Detected unwanted structural patterns"

    return True, "Identified as an epoxy fatty acid with appropriate functional groups and chain length"