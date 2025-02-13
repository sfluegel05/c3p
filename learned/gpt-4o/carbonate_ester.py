"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the structural motif of O=C(OX)OX, with variations possibly forming rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General carbonate ester pattern includes both acyclic and cyclic features
    carbonate_ester_patterns = [
        Chem.MolFromSmarts("O=C(O[*])O[*]"),  # Linear or branched with two alkoxy groups
        Chem.MolFromSmarts("O=C1OC(=O)O1"),  # Cyclic carbonate ester
        Chem.MolFromSmarts("O=C(Oc1ccccc1)Oc1ccccc1")  # Include esters with aromatic esters
    ]
    
    for pattern in carbonate_ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbonate ester structure"

    return False, "No carbonate ester pattern found"

# Example usage
# print(is_carbonate_ester("COC(=O)OC"))  # Example for dimethyl carbonate