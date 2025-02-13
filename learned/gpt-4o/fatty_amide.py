"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid, characterized by
    amide groups and typically long hydrocarbon chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary, secondary, and tertiary amide group patterns
    amide_patterns = [
        Chem.MolFromSmarts("C(=O)N"),      # Primary amide
        Chem.MolFromSmarts("C(=O)NC"),     # Secondary amide
        Chem.MolFromSmarts("C(=O)N(C)C")   # Tertiary amide
    ]

    # Check if any amide pattern matches
    amide_found = any(mol.HasSubstructMatch(pattern) for pattern in amide_patterns)
    if not amide_found:
        return False, "No amide group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count the number of double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)

    # Verify if it has characteristics of a fatty acid chain (long, possibly unsaturated)
    if c_count < 8:
        return False, f"Insufficient carbon count for fatty chain, found {c_count} carbons"
    elif c_count > 8 and double_bond_count > 0:
        return True, "Contains amide group and features of a fatty acid chain (unsaturation)"
    elif c_count >= 10:
        return True, "Contains amide group and long hydrocarbon chain characteristic of fatty acids"

    return False, "Does not match typical characteristics of fatty amides"