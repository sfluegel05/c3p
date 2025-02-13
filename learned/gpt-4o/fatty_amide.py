"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid, 
    characterized by amide groups and typically long aliphatic hydrocarbon chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary and secondary amide group patterns
    amide_pattern = Chem.MolFromSmarts("C(=O)N")

    # Check if any amide pattern matches
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False)

    # Check for a sufficiently long carbon chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCC")
    has_long_chain = mol.HasSubstructMatch(aliphatic_chain_pattern)

    if c_count >= 12 and has_long_chain:
        return True, "Contains amide group and a sufficiently long aliphatic chain characteristic of fatty acids"

    return False, "Does not match typical characteristics of fatty amides or chain is too short"