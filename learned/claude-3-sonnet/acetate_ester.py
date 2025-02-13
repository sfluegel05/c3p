"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:35676 acetate ester
Any carboxylic ester where the carboxylic acid component is acetic acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acetate_ester(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is an acetate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for acetate ester pattern: -O-C(=O)C
    acetate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX3]")
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)

    if not acetate_matches:
        return False, "No acetate ester group found"

    # Check if the carbonyl carbon is connected to a methyl group
    for match in acetate_matches:
        c_idx = match[1]  # Index of the carbonyl carbon
        c_atom = mol.GetAtomWithIdx(c_idx)
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetSymbol() == "C" and sum(n.GetSymbol() == "H" for n in neighbor.GetNeighbors()) == 3:
                return True, "Contains an acetate ester group (-O-C(=O)C)"

    return False, "No acetate ester group found"