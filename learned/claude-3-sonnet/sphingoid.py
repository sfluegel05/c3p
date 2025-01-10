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
    Sphingoids are sphinganine, its homologs and stereoisomers, and their hydroxy/unsaturated derivatives.
    Excludes N-acylated derivatives (ceramides) and complex sphingolipids.

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

    # Count basic atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Basic requirements
    if c_count < 12:
        return False, "Carbon chain too short for sphingoid"
    if n_count == 0:
        return False, "No amino group found"
    if o_count == 0:
        return False, "No hydroxyl groups found"

    # Exclude N-acylated derivatives (ceramides)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]C(=O)C")):
        return False, "N-acylated derivative (ceramide)"

    # Core structure patterns for different sphingoid types
    core_patterns = [
        # Standard sphingoid core (2-amino-1,3-diol)
        "[CH2][CH]([NH2,NH3+])[CH](O)",
        # 3-dehydro variant
        "[CH2][CH]([NH2,NH3+])C(=O)",
        # 1-deoxy variant with methyl
        "[CH3][CH]([NH2,NH3+])[CH](O)",
        # Phytosphingosine-like core
        "[CH2][CH]([NH2,NH3+])[CH](O)[CH](O)",
        # N,N-dimethyl variant
        "[CH2][CH]([NH+](C)C)[CH](O)",
        # Branched variants
        "[CH2][CH]([NH2,NH3+])[CH](O)C(C)C"
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
            
    if not has_core:
        return False, "No characteristic sphingoid core structure found"

    # Look for long carbon chain
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # At least 8 carbons in a chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain found"

    # Identify specific features
    has_double_bond = mol.HasSubstructMatch(Chem.MolFromSmarts("C=C"))
    has_extra_oh = mol.HasSubstructMatch(Chem.MolFromSmarts("[CH](O)[CH](O)"))
    has_phosphocholine = mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]OP(=O)([O-])OCC[N+](C)(C)C"))
    
    # Build classification reason
    features = []
    if has_double_bond:
        features.append("unsaturation")
    if has_extra_oh:
        features.append("additional hydroxylation")
    if has_phosphocholine:
        features.append("phosphocholine group")
        
    reason = "Contains sphingoid core structure"
    if features:
        reason += " with " + ", ".join(features)

    return True, reason