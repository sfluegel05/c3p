"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for three or more isolated double bonds
    double_bond_pattern = Chem.MolFromSmarts("[CD2]=[CD2]")  # Exclude the possibility of conjugation
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, f"Found {len(double_bond_matches)} isolated double bonds, need at least 3 for trienoic structure"
    
    # Verify that the double bonds are integrated into a long linear chain
    long_chain_pattern = Chem.MolFromSmarts("CCCCCC[CD2]=[CD2]CCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Double bonds are not part of a long carbon chain structure typical of trienoic fatty acids"
    
    # Count number of carbon atoms to ensure a long chain representing a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too short hydrocarbon chain for a typical trienoic fatty acid"

    # Check for aliphatic character and avoid excessive branching
    if not any(atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 for atom in mol.GetAtoms()):
        return False, "Molecule does not have enough sp3 hybridized carbon atoms characteristic of fatty acids"
    
    return True, "Contains carboxylic acid group and three isolated double bonds in long carbon chain typical of trienoic fatty acid"