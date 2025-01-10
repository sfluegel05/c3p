"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to the molecule if not already present
    mol = Chem.AddHs(mol)

    # Check if the molecule is acyclic
    if rdMolDescriptors.CalcNumRings(mol) != 0:
        return False, "Molecule contains rings"

    # Count the number of carbon-carbon double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and \
           bond.GetBeginAtom().GetAtomicNum() == 6 and \
           bond.GetEndAtom().GetAtomicNum() == 6:
            double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Found {double_bond_count} carbon-carbon double bonds, need exactly 1"

    # Check if the molecule is a hydrocarbon (only contains C and H)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
            return False, "Molecule contains non-carbon and non-hydrogen atoms"

    # Calculate the expected hydrogen count for an alkene (CnH2n)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    expected_h_count = 2 * c_count

    # Allow for slight deviations due to branching
    if abs(h_count - expected_h_count) > 2:
        return False, f"Hydrogen count ({h_count}) does not match expected count ({expected_h_count}) for an alkene"

    return True, "Acyclic hydrocarbon with exactly one carbon-carbon double bond"