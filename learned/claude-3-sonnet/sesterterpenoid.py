"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are derived from sesterterpenes (C25 terpenoids), though some
    carbons may be removed through modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count ({c_count}) outside typical range for sesterterpenoids (20-30)"

    # Check for rings (sesterterpenoids typically have multiple rings)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, f"Too few rings ({ring_count}) for a sesterterpenoid"

    # Check for presence of methyl groups (characteristic of terpenoids)
    methyl_pattern = Chem.MolFromSmarts("C[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 2:
        return False, "Too few methyl groups for a sesterterpenoid"

    # Check for unsaturation (double bonds are common in terpenoids)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_count == 0:
        return False, "No carbon-carbon double bonds found"

    # Calculate some molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for sesterterpenoids"

    # Check for branching (terpenoids are typically branched)
    branching_pattern = Chem.MolFromSmarts("[CH1,CH0](-[*])(-[*])(-[*])")
    branch_points = len(mol.GetSubstructMatches(branching_pattern))
    if branch_points < 2:
        return False, "Insufficient branching for a sesterterpenoid"

    # Look for typical cyclic terpenoid patterns
    decalin_pattern = Chem.MolFromSmarts("C1CCC2CCCCC2C1")
    if not mol.HasSubstructMatch(decalin_pattern) and ring_count < 3:
        return False, "Missing typical terpenoid ring structures"

    # If all checks pass, it's likely a sesterterpenoid
    return True, "Matches sesterterpenoid characteristics: C20-C30 skeleton, multiple rings, " \
                 "methyl groups, appropriate branching and unsaturation"