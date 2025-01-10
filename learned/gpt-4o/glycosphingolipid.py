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
    sphingoid_pattern = r"NC\(=O\)[C@@H]([C@@H](O)[C@H](O)CO)C"
    # This pattern expects an amide linkage with an alcohol group on the adjacent carbon

    glycosidic_linkage_pattern = r"O[C@H]([C@H](O)[C@H](O))C"
    # This pattern ensure the presence of a sugar linked via an O-glycosidic bond

    # Check for sphingoid part (a sphingoid backbone is often a primary amine with a long aliphatic chain)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(sphingoid_pattern)):
        return False, "No sphingoid or ceramide backbone identified"

    # Check for glycosidic linkage to a sugar moiety
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(glycosidic_linkage_pattern)):
        return False, "No glycosidic linkage to a sugar moiety (needed at O-1 position of the sphingoid)"

    return True, "Molecule contains a sphingoid backbone with a glycosidic linkage at the O-1 position to a sugar moiety"

# Example
example_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1OC1O[C@H](COP(O)(=O)O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)CCCCCCCCCCCCCCC"
print(is_glycosphingolipid(example_smiles))