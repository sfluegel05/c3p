"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:18244 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is a phospholipid with a sphingoid base backbone,
    an amide-linked fatty acid chain, and a phosphocholine head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[N+](C)(C)COCP(=O)([O-])")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Missing phosphocholine group"

    # Look for amide group connecting fatty acid chain
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide group for fatty acid chain"

    # Look for sphingoid base backbone (long chain with optional double bonds and substituents)
    sphingoid_pattern = Chem.MolFromSmarts("[NH1X3,NH2X2,NH0X1&r4][CX4H2][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2H]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        sphingoid_pattern = Chem.MolFromSmarts("[NH1X3,NH2X2,NH0X1&r4][CX4H2][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2H]")
        if not mol.HasSubstructMatch(sphingoid_pattern):
            return False, "Missing sphingoid base backbone"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short for sphingomyelin"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for sphingomyelin"
    if o_count < 5:
        return False, "Too few oxygens for sphingomyelin"
    if n_count < 1:
        return False, "Missing nitrogen for sphingoid base"

    return True, "Contains sphingoid base backbone with amide-linked fatty acid chain and phosphocholine head group"