"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is classified as a wax based on its SMILES string.
    Waxes typically consist of long-chain alcohols esterified to long-chain fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define ester pattern with a requirement for long chains on both sides
    ester_pattern_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_pattern_smarts)
    
    # Check for ester presence
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond present"

    # Helper function to determine if a chain has sufficient length & aliphatic nature
    def is_long_aliphatic_chain(chain):
        c_count = sum(1 for atom in chain.GetAtoms() if atom.GetAtomicNum() == 6)
        return c_count >= 15  # Modified threshold based on typical wax chain lengths

    # Iterate over all ester matches
    for ester in mol.GetSubstructMatches(ester_pattern):
        acyl_chain = Chem.PathToSubmol(mol, [ester[0]])  # Part before ester linkage
        alkoxy_chain = Chem.PathToSubmol(mol, [ester[2]])  # Part after ester linkage

        if is_long_aliphatic_chain(acyl_chain) and is_long_aliphatic_chain(alkoxy_chain):
            return True, "Contains a long-chain ester typical of waxes"

    return False, "Ester bond present, but lacking required long aliphatic chains"