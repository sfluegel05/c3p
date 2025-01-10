"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the 
    carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ester linkage pattern: C(=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Find ester group and analyze chains on either side
    matches = mol.GetSubstructMatches(ester_pattern)

    # Check for both sides of the ester
    for match in matches:
        # Get atoms on either side of the ester linkage
        carbonyl_carbon = match[0]
        oxygen = match[2]

        # Check for carbon chains extending from each side
        carbon_chain_1 = set([atom.GetIdx() for atom in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors() if atom.GetAtomicNum() == 6])
        carbon_chain_2 = set([atom.GetIdx() for atom in mol.GetAtomWithIdx(oxygen).GetNeighbors() if atom.GetAtomicNum() == 6])

        # Exclude ester group and check length of chains
        carbon_chain_1.discard(carbonyl_carbon)
        carbon_chain_1.discard(oxygen)

        carbon_chain_2.discard(carbonyl_carbon)
        carbon_chain_2.discard(oxygen)

        # Check if both chains are sufficiently long
        if len(carbon_chain_1) >= 8 and len(carbon_chain_2) >= 8:
            return True, "Molecule contains a fatty acid ester linkage with sufficient carbon chain length"
    
    return False, "Carbon chains are too short to be considered fatty acid/alcohol"