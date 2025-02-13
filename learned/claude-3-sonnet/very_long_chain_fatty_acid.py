"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:36979 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Find the longest linear carbon chain
    longest_chain_length = 0
    for chain in Chem.FindAllSubgraphsOfLengthMToN(mol, 1, 1000):
        chain_length = len(chain)
        is_linear = all(atom.GetDegree() <= 2 for atom in chain)
        if is_linear and chain_length > longest_chain_length:
            longest_chain_length = chain_length

    if longest_chain_length <= 22:
        return False, f"Longest carbon chain has {longest_chain_length} atoms, must be > 22"

    # Check for unsaturations
    unsaturated_pattern = Chem.MolFromSmarts("[CH2]=[CH]~[CH2]~[CH2]")
    unsaturated_matches = mol.GetSubstructMatches(unsaturated_pattern)
    unsaturated_lengths = [len(match) for match in unsaturated_matches]
    max_unsaturated_length = max(unsaturated_lengths, default=0)

    if max_unsaturated_length > 4:
        return False, "Too many unsaturations for fatty acid"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    allowed_atoms = set("C H O".split())
    formula_atoms = set(formula.replace("[", "").replace("]", ""))
    if not formula_atoms.issubset(allowed_atoms):
        return False, "Only C, H, and O atoms are allowed"

    return True, "Carbon chain length > 22, with carboxylic acid group"