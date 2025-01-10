"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide contains a glucose residue linked via a glycosidic bond
    to a ceramide backbone (sphingosine or a similar sphingoid base plus a fatty acid).

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

    # SMARTS pattern for beta-D-glucose moiety (more tolerant to stereochemistry)
    glucose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](CO)O1")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"
    
    # Expanded pattern for the ceramide backbone with flexibility in chain positioning
    ceramide_pattern = Chem.MolFromSmarts("N[C@@H](CO[*])C(=O)C{8,}")  # modified to accept longer chains consistently
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No adaptable ceramide-like backbone found"
    
    # Check for the sphingosine-like structure including common features
    sphingosine_pattern = Chem.MolFromSmarts("[C@H](O)C=C[C@H](O)CC")  # capturing double bonds and hydroxyls
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine-like base structure identified"

    # Refined check for a long hydrocarbon chain beyond what ceramide contains
    chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCCCCC")  # roughly 20+ carbon chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficiently long hydrocarbon chain, typical of ceramides"

    return True, "Contains a beta-D-glucose moiety with a recognizable sphingosine-like and ceramide backbone"