"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: CHEBI:XXXXX beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid contains a beta-D-glucuronic acid moiety linked via a glycosidic bond.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for beta-D-glucuronic acid with glycosidic bond
    # Matches the core structure: [O]-C1-C(=O)O with correct stereochemistry
    glucuronic_pattern = Chem.MolFromSmarts("[#8]-[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](C(=O)O)O1")
    
    if mol.HasSubstructMatch(glucuronic_pattern):
        return True, "Contains beta-D-glucuronic acid with glycosidic bond"
    return False, "No beta-D-glucuronic acid glycoside detected"