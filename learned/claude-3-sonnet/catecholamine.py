"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: catecholamine
Definition: 4-(2-Aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] and derivatives formed by substitution
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
    
    # Handle charges - get largest fragment if salt
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        mol = max(fragments, key=lambda m: m.GetNumAtoms())
    
    # Essential structure patterns
    
    # 1,2-dihydroxy benzene (catechol) with para position specified
    # [#6]:1 represents aromatic carbon
    # The numbers ensure the hydroxyls are ortho to each other
    # The last position (4) must have a carbon attached
    catechol_pattern = Chem.MolFromSmarts("[#6]:1:[#6](-[OH1,O-]):[#6](-[OH1,O-]):[#6]:[#6](-[#6]):[#6]:1")
    
    # Ethylamine chain patterns - must include connection to ring
    ethylamine_patterns = [
        # Basic ethylamine chain
        Chem.MolFromSmarts("c1c(-CCN)cc(O)c(O)c1"),
        # Allow hydroxylation and substitution
        Chem.MolFromSmarts("c1c(-CC(O)N)cc(O)c(O)c1"),
        Chem.MolFromSmarts("c1c(-C(O)CN)cc(O)c(O)c1"),
        # Allow N-substitution
        Chem.MolFromSmarts("c1c(-CCN([H,C])[H,C])cc(O)c(O)c1"),
        # Allow alpha-carbon substitution
        Chem.MolFromSmarts("c1c(-C(C)CN)cc(O)c(O)c1"),
    ]
    
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No 1,2-dihydroxybenzene (catechol) moiety found"
    
    has_ethylamine = False
    for pattern in ethylamine_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_ethylamine = True
            break
    
    if not has_ethylamine:
        return False, "No appropriate ethylamine chain found at position 4"
    
    # Additional checks
    
    # Check molecule size (catecholamines are relatively small)
    if mol.GetNumAtoms() > 50:
        return False, "Molecule too large for a typical catecholamine"
    
    # Count aromatic rings
    aromatic_rings = mol.GetSubstructMatches(Chem.MolFromSmarts("a1aaaaa1"))
    if len(aromatic_rings) > 2:  # Allow max 2 rings (e.g., for dobutamine)
        return False, "Too many aromatic rings"
    
    # Check for problematic features that shouldn't be in catecholamines
    problematic_features = [
        Chem.MolFromSmarts("C1=NC=CN1"), # Imidazole
        Chem.MolFromSmarts("C1CCCCC1"), # Cyclohexane
        Chem.MolFromSmarts("C(=O)N"), # Amide (unless part of larger necessary structure)
        Chem.MolFromSmarts("C(=O)O[H,C]"), # Ester/Acid (unless part of larger necessary structure)
        Chem.MolFromSmarts("S"), # Sulfur-containing
        Chem.MolFromSmarts("P"), # Phosphorus-containing
    ]
    
    for feature in problematic_features:
        if feature is not None and mol.HasSubstructMatch(feature):
            matches = mol.GetSubstructMatches(feature)
            # Allow if the feature is part of a known catecholamine structure
            if len(matches) > 1:  # Multiple instances are suspicious
                return False, "Contains chemical features not typical of catecholamines"
    
    return True, "Contains 1,2-dihydroxybenzene with appropriate 4-(2-aminoethyl) substitution"