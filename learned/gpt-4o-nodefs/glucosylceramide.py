"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide contains a glucose residue linked via a glycosidic bond
    to a ceramide backbone (sphingosine plus fatty acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for beta-D-glucose moiety
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](CO)O1")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"
    
    # SMARTS pattern for typical ceramide-like structure (sphingosine with a fatty acid)
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](CO[C@H]1O[C@H](CO)O[C@H](O)[C@H]1O)C")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide-like backbone found"
    
    # Check for long hydrocarbon chains, a characteristic of ceramides
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain, typical of ceramides"

    return True, "Contains a beta-D-glucose moiety and a ceramide-like backbone"