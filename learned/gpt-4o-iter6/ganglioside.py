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

    # General pattern for ceramide and sphingosine backbone
    broad_ceramide_pattern = Chem.MolFromSmarts("NC(=O)C[C@H](O)C[C@@H](O)CCCCCCCCCCCCCC")
    
    # Sialic acid substructure pattern
    broad_sialic_acid_pattern = Chem.MolFromSmarts("C(=O)C[C@H](O)[C@@H](O)C(O)CO")
    
    # Check for ceramide backbone
    if not mol.HasSubstructMatch(broad_ceramide_pattern):
        return False, "No ceramide backbone found"
    
    # Check for presence of sialic acid residue
    sialic_acid_residues = mol.GetSubstructMatches(broad_sialic_acid_pattern)
    if len(sialic_acid_residues) < 1:
        return False, "No sialic acid residues found"
    
    return True, "Contains ceramide backbone and one or more sialic acids"