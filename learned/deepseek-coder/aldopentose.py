"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:34084 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a (potential) aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 5 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected exactly 5 carbons, found {c_count}"

    # Check for potential aldehyde group (explicit or in hemiacetal form)
    # Explicit aldehyde pattern
    explicit_aldehyde = Chem.MolFromSmarts("[CX3H1](=O)")
    # Hemiacetal pattern (carbon attached to two oxygens, one single and one double bond)
    hemiacetal_pattern = Chem.MolFromSmarts("[CX4]([OX2])([OX2])")
    
    if not (mol.HasSubstructMatch(explicit_aldehyde) or mol.HasSubstructMatch(hemiacetal_pattern)):
        return False, "No potential aldehyde group found (explicit or in hemiacetal form)"

    # Check for 4 oxygen atoms (including those in hemiacetal/hemiketal forms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 5:
        return False, f"Expected exactly 5 oxygens, found {o_count}"

    # Check molecular weight (expanded range to accommodate different forms)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130 or mol_wt > 170:
        return False, f"Molecular weight {mol_wt:.2f} is outside expected range for aldopentoses"

    # Check for at least 4 hydroxyl groups (including those in hemiacetal/hemiketal forms)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 4:
        return False, f"Expected at least 4 hydroxyl groups, found {len(hydroxyl_matches)}"

    return True, "Contains exactly 5 carbons, a potential aldehyde group, and appropriate oxygen/hydroxyl groups"