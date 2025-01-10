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
    if c_count > 45:  # Allow for some modifications but exclude large molecules
        return False, f"Too many carbons ({c_count}) for a sesterterpenoid"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:  # Exclude large glycosides and complex molecules
        return False, f"Molecular weight ({mol_wt:.1f}) too high for a sesterterpenoid"
    
    # Check for characteristic sesterterpenoid patterns
    sesterterpenoid_patterns = [
        # Ophiobolin-like core
        Chem.MolFromSmarts("C1CC2CCC3C(C2)CCC3C1"),
        # Common 5-8 membered ring patterns
        Chem.MolFromSmarts("C1CCCC2CCCC12"),
        Chem.MolFromSmarts("C1CCC2CCCC2C1"),
        # Linear sesterterpenoid chain pattern
        Chem.MolFromSmarts("CCC(C)CCCC(C)CCCC(C)CCCC(C)C"),
        # Characteristic branching patterns
        Chem.MolFromSmarts("C=C(C)CCC=C(C)CCC=C(C)C")
    ]
    
    pattern_matches = sum(1 for pattern in sesterterpenoid_patterns 
                         if pattern is not None and mol.HasSubstructMatch(pattern))
    
    if pattern_matches == 0:
        return False, "No characteristic sesterterpenoid patterns found"

    # Look for steroid core to exclude steroids
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2)CCC4CCCC34C1")
    if steroid_pattern and mol.HasSubstructMatch(steroid_pattern):
        return False, "Appears to be a steroid rather than sesterterpenoid"

    # Check for methyl groups in characteristic positions
    methyl_patterns = [
        Chem.MolFromSmarts("[CH3]-[CH2,CH1,CH0]"), # Terminal methyl
        Chem.MolFromSmarts("[CH3]-C=C"),  # Methyl attached to double bond
    ]
    
    total_methyl_count = sum(len(mol.GetSubstructMatches(pattern)) 
                            for pattern in methyl_patterns if pattern is not None)
    
    if total_methyl_count < 2:
        return False, "Too few methyl groups for a sesterterpenoid"

    # Check for unsaturation (double bonds or carbonyls)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    
    unsaturation_count = (len(mol.GetSubstructMatches(double_bond_pattern)) +
                         len(mol.GetSubstructMatches(carbonyl_pattern)))
    
    if unsaturation_count == 0:
        return False, "No unsaturation found (requires double bonds or carbonyls)"

    # Count rings but allow for linear sesterterpenoids
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count == 0 and c_count < 25:
        return False, "Linear molecule with too few carbons for sesterterpenoid"
    
    # Check for characteristic branching
    branch_pattern = Chem.MolFromSmarts("[CH1,CH0](-[*])(-[*])(-[*])")
    if branch_pattern:
        branch_count = len(mol.GetSubstructMatches(branch_pattern))
        if branch_count < 2 and ring_count < 2:
            return False, "Insufficient branching for a sesterterpenoid"

    return True, "Matches sesterterpenoid characteristics with appropriate carbon count, " \
                "structural features, and branching patterns"