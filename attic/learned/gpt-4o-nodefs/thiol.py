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

    # Define SMARTS pattern for a thiol group (covering open chains, rings, etc.)
    thiol_patterns = [
        Chem.MolFromSmarts("[SX2H]"),          # standard thiol
        Chem.MolFromSmarts("[cH][SX2H]"),      # thiol adjacent to aromatic carbon
        Chem.MolFromSmarts("[C;!R]S"),         # aliphatic thiol
    ]

    # Look for thiol group in the molecule
    for pattern in thiol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a thiol group (-SH)"

    return False, "No thiol group found"


# Example usage
smiles_examples = [
    "CCSCCS", "CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCS", "SCC(CC)C",
    "CCO[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CCS(C)(=O)=O)C(=O)OC(C)C",
    "NC(=O)CCCCC(S)CCS", "SC(CS)(C)C", "Mary had a little lamb",
    "SCCCCCCCS", "NCCS"
]

# Test the function
for sml in smiles_examples:
    result, reason = is_thiol(sml)
    print(f"SMILES: {sml} => is_thiol: {result}, Reason: {reason}")