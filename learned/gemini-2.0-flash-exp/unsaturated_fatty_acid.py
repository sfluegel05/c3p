"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O or -C(=O)[O-])
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH,O-]") # handles both neutral and anionic forms.
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Does not contain a carboxylic acid group"

    # Define SMARTS patterns for C=C and C#C bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if not (mol.HasSubstructMatch(double_bond_pattern) or mol.HasSubstructMatch(triple_bond_pattern)):
        return False, "Does not contain any C=C or C#C bond"

    # Check for a long carbon chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Does not contain a sufficiently long carbon chain"

    # Check number of carbons
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 6:
        return False, "Too few carbon atoms"

    # If all checks pass
    return True, "Contains a carboxylic acid group, a long chain and at least one C=C or C#C bond"