"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    An olefinic fatty acid is characterized by having a carboxylic acid group
    and at least one C=C (carbon-carbon double bond), typically within a long aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined detection for any variant of carboxylic acid group (-C(O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Detect carbon-carbon double bond (C=C)
    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(cc_double_bond_pattern):
        return False, "No carbon-carbon double bond found"

    # Additional check for carbon chain length to ensure it's a fatty acid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8: # Assumes fatty acids typically have more than 8 carbons
        return False, "Likely not a fatty acid due to insufficient carbon atoms"

    return True, "Contains both a carboxylic acid group and at least one C=C double bond characteristic of olefinic fatty acids"