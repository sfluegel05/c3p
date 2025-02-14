"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is a fatty acid that contains one or more rings.

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

    # Check for ring structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No ring structure found"

    # Check for carboxylic acid group and a long carbon chain using SMARTS
    # look for carboxylic acid group: C(=O)O or C(=O)[O-] (as a protonated and deprotonated group)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH,O-]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
      return False, "No carboxylic acid group found"
    
    # check for long carbon chain: minimum of 4 carbons
    long_chain_pattern = Chem.MolFromSmarts('[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]')
    if not mol.HasSubstructMatch(long_chain_pattern):
      return False, "No long carbon chain found"

    # count carbons: fatty acids typically have >4 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbon atoms for a fatty acid"
   
    return True, "Contains ring and fatty acid substructures"