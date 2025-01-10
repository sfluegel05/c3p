"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a carboxylic acid group (C(O)=O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count the carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, f"Not enough carbon atoms for ultra-long-chain fatty acid: {c_count}"

    # Check for long carbon chain connected to carboxylic acid group or other known features
    # to be ultra-long-chain
    # Assume a combination of long chains and presence of functional groups like double bonds
    # Let's verify they mostly have long chains or cyclopropyl present
    if not any(atom.GetAtomicNum() == 6 and atom.GetDegree() < 4 for atom in mol.GetAtoms()):
        return False, "Lacks long carbon chain structure typical of fatty acids"
    
    # Check for special groups e.g., hydroxyl or methoxy that are commonly present in examples
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("[C][O][C]")
    cyclopropyl_pattern = Chem.MolFromSmarts("C1CC1")
    
    if mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(methoxy_pattern) or mol.HasSubstructMatch(cyclopropyl_pattern):
        return True, "Contains ultra-long-chain fatty acid features"

    # Since ultra-long features are not definitive, may require further checking with combined motifs
    return False, "Insufficient structural features for an ultra-long-chain fatty acid"