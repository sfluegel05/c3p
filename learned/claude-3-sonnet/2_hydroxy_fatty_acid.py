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
    # [CH,CH2,CH3] ensures we're matching an aliphatic carbon with 1-3 hydrogens
    # The pattern matches C(=O)(O)-C(O)- where the second C is the alpha position
    pattern = Chem.MolFromSmarts("[CH,CH2,CH3](O)C(=O)O")
    if not mol.HasSubstructMatch(pattern):
        return False, "No 2-hydroxy carboxylic acid group found"

    # Count carbons - fatty acids typically have 4+ carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Count oxygens - should have exactly 3 (COOH + OH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 3:
        return False, "Must have exactly 3 oxygen atoms (COOH + OH)"

    # Check for primarily hydrocarbon nature (C,H,O only)
    non_cho_atoms = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in [1,6,8])
    if non_cho_atoms > 0:
        return False, "Contains atoms other than C,H,O"

    # Check for single carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if len(mol.GetSubstructMatches(carboxyl_pattern)) != 1:
        return False, "Must have exactly one carboxylic acid group"

    # Check for single hydroxy group (excluding the COOH)
    hydroxy_pattern = Chem.MolFromSmarts("[OH1]")
    if len(mol.GetSubstructMatches(hydroxy_pattern)) != 1:
        return False, "Must have exactly one hydroxy group (excluding COOH)"

    # Check for linear carbon chain nature
    # Count branching points (carbons with more than 2 non-H connections)
    num_branches = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[C^3,C^4]")))
    if num_branches > 2:  # Allow some branching for methylated fatty acids
        return False, "Too many branch points for a fatty acid"

    # Success - molecule matches all criteria for 2-hydroxy fatty acid
    return True, "Matches 2-hydroxy fatty acid pattern with appropriate chain length"