"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA has a fatty acyl group with a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for coenzyme A substructure
    # This SMARTS pattern needs refinement to match more CoA fragments accurately
    # Adjust pattern based on known structure alignment with examples
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP")  # Part of CoA backbone
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A backbone found"

    # Find the longest carbon chain in the molecule
    # Use path finding to determine the longest sequence of carbon-carbon bonds
    longest_chain = 0
    for chain in rdmolops.GetLongestChain(mol):
        if len(chain) > longest_chain:
            longest_chain = len(chain)

    # Check if longest chain exceeds 22 carbons
    # Minor error in previous: be sure to analyze using atom indexes and bond types
    if longest_chain <= 22:
        return False, f"Longest carbon chain length is {longest_chain}, not greater than C22"

    return True, f"Contains CoA backbone and fatty acyl chain length is {longest_chain}, which is greater than C22"