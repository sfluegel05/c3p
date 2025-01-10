"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the hydrocarbon chain and a carboxyl group.

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

    # Check for a carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O,H]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Count the number of C=C or C#C double/triple bonds
    bond_counts = {'double': 0, 'triple': 0}
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            bond_counts['double'] += 1
        elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            bond_counts['triple'] += 1

    # Validate there is exactly one double or triple bond
    total_unsaturation = bond_counts['double'] + bond_counts['triple']
    if total_unsaturation != 1:
        return False, f"Found {total_unsaturation} unsaturations, must be exactly one"

    # Ensure the rest of the chain is a hydrocarbon (singly bonded C atoms)
    # Count carbon atoms not involved in multiple bonds
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    unsaturated_count = 2 * bond_counts['double'] + 3 * bond_counts['triple']
    
    if (carbon_count - unsaturated_count) < 4: # At least 4 apart from unsaturation
        return False, "Insufficient length of saturated chain"

    return True, "Contains a single double/triple bond in a hydrocarbon chain with a carboxyl group"