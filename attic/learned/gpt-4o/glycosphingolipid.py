"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is defined as a carbohydrate-containing derivative of 
    a sphingoid or ceramide, with the carbohydrate linked to O-1 of the sphingoid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for glycosphingolipid components
    # Refined sphingoid pattern: long carbon chain with terminal amide linkage 
    sphingoid_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](CO)C(C)O")  # More complete backbone with amide and possible hydroxyl

    # Improved sugar moiety pattern: carbohydrate ring structure with glycosidic linkage
    sugar_linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1")  # Typical hexose ring

    # Ensure valid structure detection
    if sphingoid_pattern is None or sugar_linkage_pattern is None:
        return (None, None)

    # Check for sphingoid base or ceramide (improved pattern)
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid or ceramide backbone identified, missing the typical base structure."

    # Check for glycosidic linkage to a sugar moiety (improved pattern)
    if not mol.HasSubstructMatch(sugar_linkage_pattern):
        return False, "No glycosidic linkage to a sugar moiety (needed at O-1 position of the sphingoid)."

    return True, "Molecule contains a sphingoid backbone with a glycosidic linkage at the O-1 position to a sugar moiety."

# Example test for validation
example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1OC(CO)C(O)C(O)C1O)C(O)CCCCCCCCCCCCC"
print(is_glycosphingolipid(example_smiles))