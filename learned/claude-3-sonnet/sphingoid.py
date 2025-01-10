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
    their hydroxy/unsaturated derivatives, and their N-acylated forms (ceramides).

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

    # Core sphingoid patterns
    core_patterns = [
        # Basic sphingoid backbone (2-amino-1,3-diol)
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH])[C@H,C@@H](O)CCC",
        
        # Ceramide backbone
        "[CH2X4][C@H,C@@H](NC(=O)C)[C@H,C@@H](O)CCC",
        
        # Phytosphingosine-type (additional OH)
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH])[C@H,C@@H](O)[C@H,C@@H](O)CC",
        
        # 3-dehydro variant
        "[CH2X4][C@H,C@@H]([NH2,NH3+,NH])C(=O)CCC",
        
        # N,N-dimethyl variant
        "[CH2X4][C@H,C@@H]([NH+](C)C)[C@H,C@@H](O)CCC",
        
        # 1-deoxy variant
        "[CH3][C@H,C@@H]([NH2,NH3+,NH])[C@H,C@@H](O)CCC"
    ]

    # Check for core structure
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
            
    if not has_core:
        return False, "No characteristic sphingoid core structure found"

    # Exclude common false positives
    exclude_patterns = [
        # Phosphatidylserines
        "[CH2X4][CHX4](OC(=O)C)[CH2X4]OP(=O)([O-])OC[CH]([NH3+,NH2])C(=O)[O-,OH]",
        # Other phospholipids
        "[CH2X4][CHX4](OC(=O)C)[CH2X4]OP(=O)([O-])",
        # Simple amino acids
        "[NH2][CH](C(=O)O)C"
    ]
    
    for pattern in exclude_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Structure matches phospholipid or other non-sphingoid pattern"

    # Identify specific features
    features = []
    
    # Check for ceramide
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]C(=O)C")):
        features.append("N-acylated")
    
    # Check for unsaturation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("unsaturated")
    
    # Check for additional hydroxylation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH](O)[CH](O)")):
        features.append("polyhydroxylated")
        
    # Check for glycosylation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]O[CH]1O[CH][CH][CH][CH]1")):
        features.append("glycosylated")
        
    # Check for phosphate/phosphocholine
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]OP(=O)([O-,OH])")):
        features.append("phosphorylated")

    # Build reason string
    reason = "Contains sphingoid core structure"
    if features:
        reason += " (" + ", ".join(features) + ")"

    return True, reason