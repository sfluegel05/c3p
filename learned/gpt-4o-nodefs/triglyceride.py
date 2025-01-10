"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern [CH2X4][CHX4][CH2X4]
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Ester linkage pattern
    ester_pattern = Chem.MolFromSmarts("[$([OH1]-C(=O))-!@[#6])]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"

    # Check for long carbon chains (indicative of fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, f"Missing fatty acid chains, found {len(fatty_acid_matches)}"

    # Count rotatable bonds to ensure chains are long
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be typical fatty acids"

    return True, "Contains glycerol backbone with three fatty acid chains attached via ester bonds"