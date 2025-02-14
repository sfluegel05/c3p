"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a molecule containing a D-ribose sugar and a nucleobase.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the D-ribose SMARTS pattern with relative stereochemistry. This pattern is more flexible.
    ribose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@@H]([C@H](O1)CO)O)O)")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Molecule does not contain the correct D-ribose sugar"

    # Check for a nitrogen atom attached to the ribose carbon (glycosidic bond).
    glycosidic_bond_pattern = Chem.MolFromSmarts("[NX3;H0][CX4]1[C@H]([C@@H]([C@@H]([C@H](O1)CO)O)O)")
    nucleobase_pattern = Chem.MolFromSmarts("[n;H0]1cc[cn]c1") # general pattern for pyrimidines
    nucleobase_pattern2 = Chem.MolFromSmarts("[n;H0]1c[nc][nc]1") # general pattern for purines

    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "Molecule does not contain a nucleobase attached to the ribose sugar via a glycosidic bond"
   
    if not (mol.HasSubstructMatch(nucleobase_pattern) or mol.HasSubstructMatch(nucleobase_pattern2)):
        return False, "Molecule does not contain a nucleobase attached to the ribose sugar"
    

    return True, "Molecule contains a D-ribose sugar and a nucleobase attached to it."