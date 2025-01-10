"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is composed of a glycosphingolipid with one or more sialic acids.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for ceramide backbone, including variability in chain length
    ceramide_pattern = Chem.MolFromSmarts("NC(=O)C[C@H](O)C[C@@H](O)C(C)C")
    
    # General pattern for sialic acid, allowing for variability
    sialic_acid_pattern = Chem.MolFromSmarts("C(=O)[C@@H](O)C(O)CC(O)C")
    
    # Check whether there exists a ceramide backbone
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"
    
    # Check for presence of sialic acids
    sialic_acid_residues = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_residues) < 1:
        return False, "No sialic acid residues found"
    
    return True, "Contains ceramide backbone and one or more sialic acids"