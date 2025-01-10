"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:17855 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a molecule containing a beta-D-glucuronide moiety,
    which is a glucuronic acid unit attached to another molecule via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible beta-D-glucuronide substructure pattern
    # This pattern matches the core structure of beta-D-glucuronide without being overly specific about stereochemistry
    glucuronide_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)C(=O)[O-])")
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "No beta-D-glucuronide moiety found"

    # Check if the glucuronide moiety is connected to another molecule via a glycosidic bond
    # The glycosidic bond is typically an oxygen atom connected to the anomeric carbon (C1) of the glucuronide
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)C(=O)[O-])-[OX2]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found connecting the glucuronide moiety"

    # Check for the presence of a carboxylate group (deprotonated)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for the correct number of hydroxyl groups (3) on the glucuronide ring
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 3:
        return False, "Insufficient hydroxyl groups on the glucuronide ring"

    return True, "Contains a beta-D-glucuronide moiety with a glycosidic bond and a carboxylate group"