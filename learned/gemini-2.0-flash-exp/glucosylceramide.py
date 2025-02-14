"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide contains a sphingosine or sphinganine backbone, a fatty acid attached via an amide bond and a beta-D-glucose head group.

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

    # Check for beta-D-glucose
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@H](CO)O)O[C@H]1")
    if not mol.HasSubstructMatch(glucose_pattern):
         return False, "No beta-D-glucose head group found"

    # Check for sphingosine/sphinganine backbone.
    # Look for a chain of at least 12 carbons with one alcohol and one amide.
    sphingosine_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CHX4]([OX2])[CHX4]([NX3])")
    if not mol.HasSubstructMatch(sphingosine_pattern):
         return False, "No sphingosine or sphinganine backbone found"
         
    # Check for amide bond
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    if not mol.HasSubstructMatch(amide_pattern):
         return False, "No amide bond found"
    
    # Check for long chain fatty acid
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
         return False, "No fatty acid chain found"
        
    # Check rotatable bonds and carbon count
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_rotatable < 10:
        return False, "Fatty acid chain too short"
    if c_count < 20:
        return False, "Too few carbons for glucosylceramide"
    
    return True, "Contains a beta-D-glucose head group, sphingosine/sphinganine backbone and a fatty acid via an amide bond"