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
    A flavone has a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.

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

    # Define the flavone core pattern: 2-aryl-1-benzopyran-4-one
    flavone_core_pattern = Chem.MolFromSmarts("[cH1]1[cH1][cH1][cH1][cH1][cH1]c2oc(-c3[cH1][cH1][cH1][cH1][cH1][cH1]3)cc(=O)c12")
    
    # Check if the molecule matches the flavone core pattern
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No 2-aryl-1-benzopyran-4-one core found"

    # Count the number of aromatic rings to ensure the presence of the aryl group
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, "Not enough aromatic rings for a flavone structure"

    # Check for the presence of the carbonyl group (C=O) in the benzopyran-4-one core
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found in the benzopyran-4-one core"

    return True, "Contains 2-aryl-1-benzopyran-4-one core with substituted derivatives"