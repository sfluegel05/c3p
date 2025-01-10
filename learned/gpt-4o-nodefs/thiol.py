"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol typically contains the -SH group connected to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a thiol group (must be connected to carbon)
    thiol_pattern = Chem.MolFromSmarts("[CX4,CX3;!R][SX2H]")  # upgrade aliphatic context only
    
    # Search for the thiol group in the molecule structure
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains a thiol group (-SH) attached to a carbon"

    return False, "No thiol group found"

# Example usage
smiles_examples = [
    "CCSCCS", "CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCS", "SCC(CC)C",
    "NC(=O)CCCCC(S)CCS", "SC(CS)(C)C", "NCCS", "OC(=O)C(=O)CS", 
    "SC1=CC=C(F)C=C1", "SCCCCCCCS"
]

# Test the function
for sml in smiles_examples:
    result, reason = is_thiol(sml)
    print(f"SMILES: {sml} => is_thiol: {result}, Reason: {reason}")