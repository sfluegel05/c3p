"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is characterized as an acyclic branched or unbranched hydrocarbon
    having the formula CnH2n+2, consisting entirely of hydrogen atoms and 
    saturated carbon atoms.
    
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
    
    # Check for presence of atoms other than carbon and hydrogen
    elements = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if elements.difference({6, 1}):  # Atomic number 6 is C, 1 is H
        return False, "Contains atoms other than carbon and hydrogen"
    
    # Calculate molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    
    # Expected pattern for alkanes: "C<num>C<num>H<num>H<num>"
    if not formula.startswith('C') or 'H' not in formula:
        return False, "Not a typical alkane molecular formula"
    
    # Check if the molecule is saturated (only single bonds)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            return False, "Contains unsaturated bonds (double or triple bonds present)"
    
    # Check for acyclic structure (cannot have rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, not acyclic"
    
    # Check if the formula matches CnH2n+2
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if hydrogen_count != 2 * carbon_count + 2:
        return False, f"Formula C{carbon_count}H{hydrogen_count} does not match CnH2n+2"

    return True, "Molecule matches the definition of an alkane"