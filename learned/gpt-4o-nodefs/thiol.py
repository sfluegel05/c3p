"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol typically contains the -SH group connected to a carbon, 
    but can also be attached as part of different molecular contexts.

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

    # Define a more comprehensive SMARTS pattern for a thiol group
    thiol_patterns = [
        Chem.MolFromSmarts("[CX4,CX3;!R][SX2H]"),  # Aliphatic thiol
        Chem.MolFromSmarts("[c][SX2H]"),           # Aromatic thiol
        Chem.MolFromSmarts("[SX2H]")               # General -SH group
    ]
    
    # Check for the presence of any thiol pattern
    for pattern in thiol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a thiol group (-SH)"

    return False, "No thiol group found"

# Example usage (for controlled testing)
smiles_examples = [
    "CCSCCS", "CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCS", "SCC(CC)C",
    "NC(=O)CCCCC(S)CCS", "SC(CS)(C)C", "NCCS", "OC(=O)C(=O)CS", 
    "SC1=CC=C(F)C=C1", "SCCCCCCCS"
]

# Test the function
for sml in smiles_examples:
    result, reason = is_thiol(sml)
    print(f"SMILES: {sml} => is_thiol: {result}, Reason: {reason}")