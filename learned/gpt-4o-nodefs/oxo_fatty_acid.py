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
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for additional ketone group (not part of carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[C;!$(C(=O)O)]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No additional ketone group found"

    # Improved check for a sufficient carbon skeleton
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # more specific for a longer alkane chain
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon backbone"
    
    # Unsaturation checking with multiple bonds or rings
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(unsaturation_pattern):
        cyclic_pattern = Chem.MolFromSmarts("C1CCCCC1")  # check if cyclic when double bonds aren't present
        if not mol.HasSubstructMatch(cyclic_pattern):
            return False, "Lacks unsaturation or cyclic components typical of oxo fatty acids"

    # Considering stereochemistry by tracking stereocenters
    if any(atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()):
        return True, "Contains carboxylic acid, ketone group, and a carbon skeleton with stereochemistry"

    return True, "Contains carboxylic acid, additional ketone group, and a sufficient carbon backbone"