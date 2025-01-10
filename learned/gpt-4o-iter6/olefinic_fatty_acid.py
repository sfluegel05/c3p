"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    An olefinic fatty acid is characterized by having a carboxylic acid group
    and at least one C=C (carbon-carbon double bond) in an aliphatic chain,
    with certain carbon chain characteristics typical of fatty acids.

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
    
    # Look for carboxylic acid group - allow for common variations
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")  # Include variations like -COO- and -COOH
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group or valid variant found"

    # Look for C=C double bonds in an aliphatic environment
    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(cc_double_bond_pattern):
        return False, "No non-ring carbon-carbon double bond found"
    
    # Verify the linearity and length of the chain
    # Ensure there is a continuous chain containing both the carboxylic group and a C=C bond.
    carbon_chain_pattern = Chem.MolFromSmarts("C=CCCCCCCCC(=O)O")  # Basic long chain with C=C and COOH
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No sufficiently extended carbon chain with both features found"

    # Count total number of carbon atoms to verify typical fatty acid length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:  # Raising the bar for typical fatty acids' length
        return False, "Carbon count too low for typical fatty acids"
    
    # Check for excessive functionalization and non-typical structures
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in (7, 8, 9, 15, 16))
    if heteroatom_count > 5:
        return False, "Excessive heteroatoms for a typical fatty acid"

    return True, "Contains both a carboxylic acid group and at least one aliphatic C=C double bond typical of olefinic fatty acids"