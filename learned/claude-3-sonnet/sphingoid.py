"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids include sphinganine, its homologs and stereoisomers, 
    their hydroxy/unsaturated derivatives, and their N-acylated forms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core sphingoid patterns with more specific requirements
    core_patterns = [
        # Basic sphingoid backbone (2-amino-1,3-diol) with various stereochemistry
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH,N])[C@H,C@@H](O)CCCCCCC",
        
        # Sphingosine-type with double bond
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH,N])[C@H,C@@H](O)CC\\C=C\\C",
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH,N])[C@H,C@@H](O)CC/C=C/C",
        
        # N-methylated variants
        "[CH2X4][C@H,C@@H]([NH+](C)C)[C@H,C@@H](O)CCCCCCC",
        
        # 3-dehydro variant
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH,N])C(=O)CCCCCCC"
    ]

    # Check for core structure
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
            
    if not has_core:
        return False, "No sphingoid core structure found"

    # Count carbons and check chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:  # Most sphingoids have 14+ carbons
        return False, "Carbon chain too short for sphingoid"

    # Exclude false positives
    exclude_patterns = [
        # Complex branched glycosides
        "[CH2]O[CH]1O[CH][CH][CH][CH]1O[CH]2O[CH][CH][CH][CH]2",
        # Peptides
        "[NH2][CH]([CH2,CH3])C(=O)[NH][CH]([CH2,CH3])C(=O)",
        # Cyclic structures with multiple nitrogens
        "C1[NH]C[NH]C1",
        # Complex oligosaccharides
        "[CH]1O[CH][CH]O[CH][CH]O[CH]O1"
    ]
    
    for pattern in exclude_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Structure matches non-sphingoid pattern"

    # Check for characteristic features
    features = []
    
    # Check for unsaturation (double bonds)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("unsaturated")
    
    # Check for N-acylation (ceramides)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]C(=O)C")):
        features.append("N-acylated")
    
    # Check for additional hydroxylation
    hydroxy_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    if hydroxy_count > 2:
        features.append("polyhydroxylated")
    
    # Check for simple glycosylation (not complex oligosaccharides)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]O[CH]1O[CH]([CH]([CH]([CH]1O)O)O)CO")):
        features.append("glycosylated")
        
    # Check for phosphorylation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]OP(=O)([O-,OH])")):
        features.append("phosphorylated")

    # Validate overall composition
    if c_count > 50:  # Too large for typical sphingoid
        return False, "Molecule too large for sphingoid"

    # Build reason string
    reason = "Contains sphingoid core structure"
    if features:
        reason += " (" + ", ".join(features) + ")"

    return True, reason