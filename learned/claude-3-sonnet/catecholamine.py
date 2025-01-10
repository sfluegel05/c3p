"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: catecholamine
Definition: 4-(2-Aminoethyl)pyrocatechol and derivatives formed by substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Convert to neutral form if possible (handle salt forms)
    mol = Chem.RemoveHs(mol)  # Remove explicit hydrogens
    
    # More flexible catechol patterns allowing different substitutions
    catechol_patterns = [
        Chem.MolFromSmarts("c1(O)c(O)cccc1"),  # 1,2-dihydroxy
        Chem.MolFromSmarts("c1c(O)c(O)ccc1"),   # 2,3-dihydroxy
        Chem.MolFromSmarts("c1cc(O)c(O)cc1"),   # 3,4-dihydroxy
        Chem.MolFromSmarts("c1c(O)cc(O)cc1"),   # 2,4-dihydroxy
        Chem.MolFromSmarts("c1(O)c(O)c(O)ccc1"), # Trihydroxy variants
        Chem.MolFromSmarts("c1(O)c(O)cc(O)cc1"),
    ]
    
    has_catechol = False
    for pat in catechol_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_catechol = True
            break
            
    if not has_catechol:
        return False, "No catechol moiety found"
    
    # More comprehensive ethylamine patterns
    ethylamine_patterns = [
        # Basic patterns
        Chem.MolFromSmarts("[cR1]!@CC[NH2,NH1,NH0]"), # Direct ethylamine chain
        Chem.MolFromSmarts("[cR1]!@CC([*,H])N"),      # Allow substitution at beta carbon
        Chem.MolFromSmarts("[cR1]!@C([*,H])CN"),      # Allow substitution at alpha carbon
        # Hydroxylated variants
        Chem.MolFromSmarts("[cR1]!@CC(O)N"),          # Beta-hydroxyl
        Chem.MolFromSmarts("[cR1]!@C(O)CN"),          # Alpha-hydroxyl
        # Various amine substitutions
        Chem.MolFromSmarts("[cR1]!@CCN([*,H])[*,H]"), # Secondary amine
        Chem.MolFromSmarts("[cR1]!@CCN([*,H])([*,H])"), # Tertiary amine
    ]
    
    has_ethylamine = False
    for pat in ethylamine_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_ethylamine = True
            break
            
    if not has_ethylamine:
        return False, "No ethylamine chain found"
    
    # Additional checks to avoid false positives
    
    # Check molecule size (catecholamines are relatively small)
    if mol.GetNumAtoms() > 40:  # Increased from 30 to allow for some larger derivatives
        return False, "Molecule too large for a typical catecholamine"
    
    # Check for maximum one benzene ring with catechol pattern
    aromatic_rings = mol.GetSubstructMatches(Chem.MolFromSmarts("a1aaaaa1"))
    if len(aromatic_rings) > 2:  # Allow max 2 rings for structures like dobutamine
        return False, "Too many aromatic rings for a catecholamine"
        
    # Check for presence of unlikely ring systems in catecholamines
    complex_ring_systems = [
        Chem.MolFromSmarts("C1=CC2=C(C=C1)C1=C(C=CC=C1)N2"), # Complex fused rings
        Chem.MolFromSmarts("C1CCCCC1"), # Cyclohexane
        Chem.MolFromSmarts("C1CCCC1"),  # Cyclopentane
    ]
    
    for ring in complex_ring_systems:
        if ring is not None and mol.HasSubstructMatch(ring):
            return False, "Contains ring systems not typical for catecholamines"
    
    return True, "Contains catechol moiety with ethylamine chain in appropriate arrangement"