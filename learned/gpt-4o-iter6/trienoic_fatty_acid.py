"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds.

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
    
    # Check for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Find and count all C=C double bonds
    double_bond_matches = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
                           bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6]
    if len(double_bond_matches) != 3:
        return False, f"Found {len(double_bond_matches)} carbon-carbon double bonds, need exactly 3"
    
    # Ensure double bonds are part of a significant carbon chain
    # Assume that "long carbon chain" means a span of at least ~12 carbons including double bonds
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if num_carbons < 12:
        return False, f"Insufficient number of carbon atoms ({num_carbons}) for a standard trienoic fatty acid"

    # Optional: Examine stereochemistry if relevant to the definition of the fatty acid
    # This could check for the specific geometric configuration using properties of the double bonds

    return True, "Contains carboxylic acid group and three double bonds in a sufficient carbon chain"