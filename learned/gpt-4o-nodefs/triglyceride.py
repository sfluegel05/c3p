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
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone pattern [C(CO)C(O)C(O)]
    glycerol_pattern = Chem.MolFromSmarts("C(CO)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for ester groups [C(=O)O]
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"
    
    # Check for long carbon chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("C~C~C~C~C")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, f"Missing fatty acid chains, found {len(fatty_acid_matches)}"
    
    # Count the number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains not long enough to represent typical fatty acids"
    
    return True, "Contains glycerol backbone with three fatty acid chains connected via ester bonds"