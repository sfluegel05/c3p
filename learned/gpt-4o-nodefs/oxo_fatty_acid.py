"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for an additional ketone group (not part of carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("[CX3]=O")
    matches = mol.GetSubstructMatches(ketone_pattern)
    if len(matches) < 2:
        return False, "No additional ketone group found"

    # Refined check for sufficient carbon backbone
    # Require long chain of carbons (typical for fatty acids), optionally with unsaturation
    carbon_chain_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient or interrupted carbon backbone"

    # Check for unsaturation or cyclic structures
    unsaturation_pattern = Chem.MolFromSmarts("C=C")  # double bond for unsaturation
    cyclic_pattern = Chem.MolFromSmarts("C1CCCCC1")  # typical 6-membered ring
    if not mol.HasSubstructMatch(unsaturation_pattern) and not mol.HasSubstructMatch(cyclic_pattern):
        return False, "Lacks unsaturation or cyclic components typical of oxo fatty acids"

    # Check stereochemistry if relevant
    chiral_centers = [atom.GetChiralTag() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED]
    if chiral_centers:
        return True, "Contains carboxylic acid, ketone group, and a carbon skeleton with stereochemistry"

    return True, "Contains carboxylic acid, additional ketone group, and a sufficient carbon backbone"