"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain fatty acid with exactly 18 carbons, two C=C double bonds,
    and a carboxylic group at one end.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 18 carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 18:
        return False, f"Has {carbon_count} carbon atoms, expected 18"

    # Look for two C=C double bonds along a continuous chain
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Count noble atoms as bonds connecting C atoms
    c_bond_count = sum(1 for atom in c_atoms if atom.GetDegree() >= 3)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    
    if double_bond_count != 2 or c_bond_count > 0:
        return False, f"Has {double_bond_count} C=C double bonds on a non-linear chain; expected exactly 2 on a linear chain"

    # Check for carboxylic acid group, which is an oxygen double bonded to a carbon (C=O) and a hydroxyl group (O-H) on the same carbon
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    # If all checks pass, it's an octadecadienoic acid
    return True, "Contains an 18-carbon straight chain, 2 C=C double bonds, and a carboxylic acid group"

# Sample usage
smiles = "CCCCCC\C=C\C=C/CCCCCCCC(O)=O"
print(is_octadecadienoic_acid(smiles))