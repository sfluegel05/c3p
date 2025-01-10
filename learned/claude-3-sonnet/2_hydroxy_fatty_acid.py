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

    # Look for 2-hydroxy carboxylic acid pattern
    # Match both protonated and deprotonated forms
    # The C-O bond must be single bond (not part of C=O)
    pattern = Chem.MolFromSmarts("[$([CX4;!r][OX2H1]),$([CX4;!r][OX1-])]C([OX2H1,OX1-])=O")
    if not mol.HasSubstructMatch(pattern):
        return False, "No 2-hydroxy carboxylic acid group found"

    # Check for aromatic atoms - fatty acids should not have any
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms - not a fatty acid"

    # Count carbons - fatty acids typically have 4+ carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short to be a fatty acid"

    # Check for primarily hydrocarbon nature (C,H,O only)
    non_cho_atoms = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in [1,6,8])
    if non_cho_atoms > 0:
        return False, "Contains atoms other than C,H,O"

    # Check for single carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1-]")
    if len(mol.GetSubstructMatches(carboxyl_pattern)) != 1:
        return False, "Must have exactly one carboxylic acid group"

    # Verify the OH is at alpha position
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[CX4;!r]([OX2H1])C([OX2H1,OX1-])=O")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "Hydroxy group not at alpha position"

    # Check carbon chain characteristics
    # Allow for branching but require a main chain of at least 4 carbons
    chain_pattern = Chem.MolFromSmarts("C~C~C~C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient carbon chain found"

    # Success - molecule matches all criteria for 2-hydroxy fatty acid
    return True, "Matches 2-hydroxy fatty acid pattern with appropriate chain length"