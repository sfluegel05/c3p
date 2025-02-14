"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid that contains three double bonds.

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

    # Check for terminal carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH1]")  # Ensure it's a standard carboxylic acid
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # Find distinct C=C double bonds along the main chain
    double_bond_pattern = Chem.MolFromSmarts("C=CC")  # Pattern including adjacent carbons
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    unique_double_bonds = set()
    for match in double_bond_matches:
        # Ensure bonds are part of the main chain, not branches
        unique_double_bonds.add(tuple(sorted((match[0], match[1]))))

    if len(unique_double_bonds) != 3:
        return False, f"Found {len(unique_double_bonds)} distinct double bonds along chain, need exactly 3"
    
    # Check for acyclic structure in the main chain context
    if Chem.rdchem.Mol(rdmol=mol).GetRingInfo().IsAromatic():
        return False, "Fatty acids should be acyclic, likely aromatic features detected"

    # Ensure sufficient length for a fatty acid
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() <= 3)
    if carbon_chain_length < 12:
        return False, "Carbon chain too short for a typical fatty acid"

    return True, "Contains an acyclic fatty acid structure with three distinct double bonds"