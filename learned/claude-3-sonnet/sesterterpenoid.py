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
    if c_count > 70:  # Increased limit to allow for complex derivatives
        return False, f"Too many carbons ({c_count}) for a sesterterpenoid"

    # Calculate molecular weight - increased limit
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1500:  # Increased to allow for glycosylated derivatives
        return False, f"Molecular weight ({mol_wt:.1f}) too high for a sesterterpenoid"
    
    # Enhanced sesterterpenoid patterns
    sesterterpenoid_patterns = [
        # Basic ophiobolin-like cores
        Chem.MolFromSmarts("C1CC2CCC3C(C2)CCC3C1"),
        Chem.MolFromSmarts("C1CC2CCC3(C)C(C2)CCC3C1"),
        # Scalarin-like patterns
        Chem.MolFromSmarts("C1CC2CCC3(C)C(C2)CC(=O)OC3C1"),
        # Common ring patterns
        Chem.MolFromSmarts("C1CCCC2CCCC12"),
        Chem.MolFromSmarts("C1CCC2CCCC2C1"),
        Chem.MolFromSmarts("C1CC2CCCC2CC1"),
        # Linear and branched patterns
        Chem.MolFromSmarts("CCC(C)CCCC(C)CCCC(C)C"),
        Chem.MolFromSmarts("C=C(C)CCC=C(C)CCC=C(C)C"),
        # Characteristic side chains
        Chem.MolFromSmarts("CC(C)=CCCC(C)=CCCC(C)=C"),
        # Ophiobolin-specific patterns
        Chem.MolFromSmarts("C1CC2C(CC1)C(=O)C=C2C"),
        # Additional ring systems
        Chem.MolFromSmarts("C1CC2C3CCC(C3)C2C1"),
        Chem.MolFromSmarts("C1CC2C(C1)C1CCC2C1")
    ]
    
    pattern_matches = sum(1 for pattern in sesterterpenoid_patterns 
                         if pattern is not None and mol.HasSubstructMatch(pattern))
    
    if pattern_matches == 0:
        # Check for simpler patterns for linear sesterterpenoids
        linear_patterns = [
            Chem.MolFromSmarts("CC(C)=CCC=C(C)CCC=C(C)C"),
            Chem.MolFromSmarts("CC(=O)C=CC(C)=CCC=C(C)C"),
            Chem.MolFromSmarts("CC(O)=CCC=C(C)CCC=C(C)C")
        ]
        linear_matches = sum(1 for pattern in linear_patterns 
                           if pattern is not None and mol.HasSubstructMatch(pattern))
        if linear_matches == 0:
            return False, "No characteristic sesterterpenoid patterns found"

    # Look for steroid core to exclude steroids
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2)CCC4CCCC34C1")
    if steroid_pattern and mol.HasSubstructMatch(steroid_pattern) and pattern_matches < 2:
        return False, "Appears to be a steroid rather than sesterterpenoid"

    # Check for methyl groups
    methyl_patterns = [
        Chem.MolFromSmarts("[CH3]-[CH2,CH1,CH0]"),
        Chem.MolFromSmarts("[CH3]-C=C"),
        Chem.MolFromSmarts("[CH3]-C(-[*])(-[*])-[*]")  # Branched methyl
    ]
    
    total_methyl_count = sum(len(mol.GetSubstructMatches(pattern)) 
                            for pattern in methyl_patterns if pattern is not None)
    
    if total_methyl_count < 2:
        return False, "Too few methyl groups for a sesterterpenoid"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    
    unsaturation_count = (len(mol.GetSubstructMatches(double_bond_pattern)) +
                         len(mol.GetSubstructMatches(carbonyl_pattern)))
    
    if unsaturation_count == 0:
        return False, "No unsaturation found (requires double bonds or carbonyls)"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count == 0 and c_count < 23:  # Slightly relaxed for linear sesterterpenoids
        return False, "Linear molecule with too few carbons for sesterterpenoid"
    
    # Check for characteristic branching
    branch_pattern = Chem.MolFromSmarts("[CH1,CH0](-[*])(-[*])(-[*])")
    if branch_pattern:
        branch_count = len(mol.GetSubstructMatches(branch_pattern))
        if branch_count < 2 and ring_count < 2:
            return False, "Insufficient branching for a sesterterpenoid"

    return True, "Matches sesterterpenoid characteristics with appropriate carbon count, " \
                "structural features, and branching patterns"