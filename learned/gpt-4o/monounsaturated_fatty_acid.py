"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have exactly one double or triple bond in a long linear carbon chain with an end carboxyl group.

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

    # Check for a carboxyl group connected directly to a carbon chain
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Count the number of C=C double and C#C triple bonds for unsaturations
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    # Find and count the unsaturated bonds
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    triple_bonds = mol.GetSubstructMatches(triple_bond_pattern)
    total_unsaturations = len(double_bonds) + len(triple_bonds)

    # Validate exactly one unsaturation
    if total_unsaturations != 1:
        return False, f"Incorrect number of unsaturations: {total_unsaturations}, need exactly one"

    # Ensure the molecule is a long carbon chain ending with a carboxyl group
    # Assume a minimum chain length by ensuring more than 6 carbons in the chain
    carbon_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_atom_count < 7:
        return False, f"Insufficient number of carbon atoms for MUFA, found {carbon_atom_count}"

    return True, "Single unsaturation in a mostly linear hydrocarbon chain with a carboxyl group present"