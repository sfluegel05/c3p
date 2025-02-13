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
    Sesterterpenoids are derived from sesterterpenes (C25 terpenoids), though they
    can have additional carbons from modifications or fewer from degradation.

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
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) for a sesterterpenoid (minimum 20)"
        
    # Check for rings (sesterterpenoids typically have multiple rings)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, f"Too few rings ({ring_count}) for a sesterterpenoid"

    # Look for methyl groups (characteristic of terpenoids)
    # More comprehensive pattern that includes different types of methyl groups
    methyl_patterns = [
        Chem.MolFromSmarts("[CH3]-[CH2,CH1,CH0]"), # Terminal methyl
        Chem.MolFromSmarts("[CH3]-C=C"),  # Methyl attached to double bond
        Chem.MolFromSmarts("[CH3]-C(-[*])(-[*])") # Branched methyl
    ]
    
    total_methyl_count = sum(len(mol.GetSubstructMatches(pattern)) 
                            for pattern in methyl_patterns if pattern is not None)
    
    if total_methyl_count < 2:
        return False, "Too few methyl groups for a sesterterpenoid"

    # Check for polycyclic structures common in terpenoids
    polycyclic_patterns = [
        Chem.MolFromSmarts("C1CCC2CCCCC2C1"), # Decalin
        Chem.MolFromSmarts("C1CC2CCC1CC2"),   # Bicyclo[3.3.0]octane
        Chem.MolFromSmarts("C1CCC2CC2CC1"),   # Bicyclo[4.2.0]octane
    ]
    
    has_polycyclic = any(mol.HasSubstructMatch(pattern) 
                        for pattern in polycyclic_patterns if pattern is not None)
    
    if not has_polycyclic and ring_count < 3:
        return False, "Missing characteristic terpenoid ring structures"

    # Check for unsaturation (double bonds are common in terpenoids)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Also check for carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    carbonyl_count = len(mol.GetSubstructMatches(carbonyl_pattern))
    
    if double_bond_count + carbonyl_count == 0:
        return False, "No unsaturation found (requires double bonds or carbonyls)"

    # Check for branching (terpenoids are typically branched)
    branching_patterns = [
        Chem.MolFromSmarts("[CH1,CH0](-[*])(-[*])(-[*])"), # Tertiary/quaternary carbon
        Chem.MolFromSmarts("C=C(C)C")  # Branched alkene
    ]
    
    total_branch_points = sum(len(mol.GetSubstructMatches(pattern)) 
                            for pattern in branching_patterns if pattern is not None)
    
    if total_branch_points < 2:
        return False, "Insufficient branching for a sesterterpenoid"

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for a sesterterpenoid"

    return True, "Matches sesterterpenoid characteristics: polycyclic structure, " \
                "methyl groups, appropriate branching and unsaturation"