"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: CHEBI:36343 Cyclic fatty acid
A cyclic fatty acid is defined as any fatty acid containing anywhere in its structure a ring of atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for rings
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No rings found in the structure"

    # Check for long carbon chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Carbon chains too short to be a fatty acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 10:
        return False, "Too few carbons for a fatty acid"
    if o_count < 2:
        return False, "Too few oxygens for a fatty acid"

    return True, "Contains a carboxylic acid group and at least one ring in the structure"