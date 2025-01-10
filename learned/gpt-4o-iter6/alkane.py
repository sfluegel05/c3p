"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is characterized by the formula CnH2n+2 and consists of
    saturated carbon atoms and hydrogen atoms, without any rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of only carbon and hydrogen atoms
    elements = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if elements.difference({6, 1}):  # Atomic number 6 is C, 1 is H
        return False, "Contains atoms other than carbon and hydrogen"
    
    # Check if the molecule is saturated (only single bonds)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            return False, "Contains unsaturated bonds (double or triple bonds present)"
    
    # Check for acyclic structure (cannot have rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, not acyclic"
    
    # Correctly count carbons and hydrogens
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

    # Calculate expected hydrogen count based on alkanes' CnH2n+2 rule
    expected_hydrogen_count = 2 * carbon_count + 2

    if hydrogen_count != expected_hydrogen_count:
        # Debug potential miscounting with a reasoned explanation
        return False, f"Formula C{carbon_count}H{hydrogen_count} does not match CnH2n+2 (expected H{expected_hydrogen_count})"
    
    return True, "Molecule matches the definition of an alkane"