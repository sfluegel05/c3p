"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy group at the alpha position (2-position)
    relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Look for 2-hydroxy pattern - more general pattern that matches alpha-hydroxy acids
    # [C,c] matches any carbon (aliphatic or aromatic)
    # [OH1,OH0] matches both protonated and deprotonated hydroxy groups
    # The pattern looks for C(=O)O-C(-[OH])-C sequence
    hydroxy_acid_pattern = Chem.MolFromSmarts("[C,c](=O)[OH1,OH0]-[C,c](-[OH1,OH0])")
    if not mol.HasSubstructMatch(hydroxy_acid_pattern):
        return False, "No hydroxy group at 2-position"

    # Count carbons to ensure it's a fatty acid (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short to be a fatty acid"

    # Count oxygens - should have at least 3 (COOH + OH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms"

    # Check for primarily hydrocarbon nature
    non_cho_atoms = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in [1,6,8])
    if non_cho_atoms > 0:
        return False, "Contains unexpected atoms"

    # Additional check to ensure the hydroxy group is actually at position 2
    # This uses a more specific pattern that ensures the OH is on the carbon
    # adjacent to the carboxylic acid
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[CH1,CH2](-[OH1])(-[C,c](=O)[OH1,OH0])")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "Hydroxy group not at alpha position"

    # Success - found all required patterns
    return True, "Contains carboxylic acid with hydroxy group at 2-position and appropriate carbon chain"