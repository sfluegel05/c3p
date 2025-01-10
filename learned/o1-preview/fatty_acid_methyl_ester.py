"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:XXXX fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester obtained by the formal condensation
    of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl ester SMARTS pattern: carbonyl carbon double bonded to oxygen,
    # single bonded to an oxygen connected to a methyl group
    methyl_ester_smarts = "[CX3](=O)[O][CH3]"
    methyl_ester_pattern = Chem.MolFromSmarts(methyl_ester_smarts)

    # Find all methyl ester groups in the molecule
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    num_methyl_esters = len(methyl_ester_matches)

    if num_methyl_esters != 1:
        return False, f"Found {num_methyl_esters} methyl ester groups, need exactly 1"

    # Define the ester SMARTS pattern: carbonyl carbon double bonded to oxygen,
    # single bonded to an oxygen connected to any carbon
    ester_smarts = "[CX3](=O)[O][#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)

    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_esters = len(ester_matches)

    if num_esters != 1:
        return False, f"Found {num_esters} ester groups, need exactly 1"

    return True, "Molecule is a fatty acid methyl ester"