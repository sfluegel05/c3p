"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is a sphingoid base with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the sphingoid base pattern (long-chain amino alcohol)
    sphingoid_base_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][OX2H]")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    # Look for the amide-linked fatty acid pattern
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) == 0:
        return False, "No amide-linked fatty acid found"

    # Check the chain length of the fatty acid (14 to 26 carbon atoms)
    fatty_acid_chain = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain)
    if len(fatty_acid_matches) < 1:
        return False, "Fatty acid chain too short"

    # Count the number of carbon atoms in the fatty acid chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 26:
        return False, f"Fatty acid chain length {c_count} is not within 14 to 26 carbon atoms"

    # Check for the presence of a hydroxyl group on carbon 2
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][CX4][OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group on carbon 2 found"

    return True, "Contains sphingoid base with amide-linked fatty acid and hydroxyl group on carbon 2"