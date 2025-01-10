"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:72010 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible flavone core pattern: 2-aryl-1-benzopyran-4-one
    # This pattern allows for substitutions on the aryl group and the benzopyran-4-one core
    flavone_core_pattern = Chem.MolFromSmarts("[c;H1]1[c;H1][c;H1][c;H1][c;H1][c;H1]1-[c;H1]1[c;H1]2[c;H1]([o;H0][c;H0](=O)[c;H1][c;H1]2)[c;H1][c;H1][c;H1]1")
    
    # Check if the molecule contains the flavone core
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No 2-aryl-1-benzopyran-4-one core found"

    # Check for the presence of the carbonyl group (C=O) in the benzopyran-4-one part
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 1:
        return False, "No carbonyl group found in the benzopyran-4-one core"

    # Check for the presence of the aryl group (benzene ring) at the 2-position
    aryl_pattern = Chem.MolFromSmarts("c1ccccc1")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern)
    if len(aryl_matches) < 1:
        return False, "No aryl group found at the 2-position"

    # Check for the presence of the oxygen in the benzopyran ring
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 1:
        return False, "No oxygen found in the benzopyran ring"

    # If all checks pass, the molecule is a flavone
    return True, "Contains 2-aryl-1-benzopyran-4-one core structure"