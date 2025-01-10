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
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Look for 2-hydroxy pattern (OH on carbon adjacent to COOH)
    hydroxy_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]-[CX4;H1,H2]([OX2H1])")
    if not mol.HasSubstructMatch(hydroxy_acid_pattern):
        return False, "No hydroxy group at 2-position"

    # Count carbons to ensure it's a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short to be a fatty acid"

    # Count oxygens - should have at least 3 (COOH + OH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms"

    # Additional check for aliphatic chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No aliphatic chain found"

    # Check that molecule isn't too complex (should be primarily hydrocarbon chain)
    other_atoms = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() not in [1,6,8])
    if other_atoms > 0:
        return False, "Contains unexpected atoms"

    # Success - found all required patterns
    return True, "Contains carboxylic acid with hydroxy group at 2-position and appropriate carbon chain"