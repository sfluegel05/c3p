"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the long carbon chain and a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a MUFA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a carboxyl group (acetic acid functional group)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Count the number of C=C double and C#C triple bonds for unsaturations
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    triple_bonds = mol.GetSubstructMatches(triple_bond_pattern)
    total_unsaturations = len(double_bonds) + len(triple_bonds)

    # Validate exactly one double or triple bond
    if total_unsaturations != 1:
        return False, f"Incorrect number of unsaturations: {total_unsaturations}, need exactly one"

    # Check for a sufficient length of hydrocarbon chain - assume minimum chain length of 8 for MUFAs
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 8:
        return False, f"Insufficient carbon chain length: {len(carbon_atoms)} atoms, minimum is 8"

    return True, "Single unsaturation in a long hydrocarbon chain with a carboxyl group present"