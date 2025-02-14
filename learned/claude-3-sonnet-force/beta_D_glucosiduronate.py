"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:36937 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a carbohydrate acid derivative anion obtained by deprotonation
    of the carboxy group of any beta-D-glucosiduronic acid.

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
    
    # Look for glucuronide substructure (glucose with carboxylate group)
    glucuronide_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)C([O-])=O")
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "No glucuronide substructure found"
    
    # Check for minimum number of hydroxy groups (at least 3)
    num_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if num_hydroxy < 3:
        return False, f"Expected at least 3 hydroxy groups, got {num_hydroxy}"
    
    # Ensure glucuronide substructure is part of a larger glucose-like structure
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)[C@@H]([C@H]([C@@H]([C@H](O)O)O)O)O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "Glucuronide substructure not part of a larger glucose-like structure"
    
    return True, "Contains glucuronide substructure (glucose with carboxylate group)"