"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with exactly one carbon-carbon double bond and formula CnH2n.

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

    # Check for hydrocarbon (only C and H atoms)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {1, 6}:
            return False, "Contains non-carbon/hydrogen atoms"

    # Count double bonds and check for triple bonds
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds += 1
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            return False, "Contains triple bonds"

    if double_bonds != 1:
        return False, f"Found {double_bonds} double bonds, need exactly 1"

    # Check for acyclic structure (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"

    # Add hydrogens to get accurate hydrogen count
    mol_with_h = Chem.AddHs(mol)
    
    # Verify molecular formula CnH2n (including all hydrogens)
    c_count = sum(1 for atom in mol_with_h.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol_with_h.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if h_count != 2 * c_count:
        return False, f"Molecular formula C{c_count}H{h_count} does not match CnH2n"

    return True, "Acyclic hydrocarbon with exactly one carbon-carbon double bond and formula CnH2n"