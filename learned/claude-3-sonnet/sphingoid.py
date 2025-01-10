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
    Sphingoids are characterized by a long carbon chain with an amino group
    and at least one hydroxyl group.

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
    if c_count < 12:  # Need a long carbon chain
        return False, "Carbon chain too short for sphingoid"
    if n_count == 0:  # Must have at least one nitrogen
        return False, "No amino group found"
    if o_count == 0:  # Must have at least one oxygen
        return False, "No hydroxyl groups found"

    # Look for primary alcohol (-CH2OH) or phosphocholine
    alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    phosphocholine = Chem.MolFromSmarts("[CH2]OP(=O)([O-])OCC[N+](C)(C)C")
    if not (mol.HasSubstructMatch(alcohol_pattern) or mol.HasSubstructMatch(phosphocholine)):
        return False, "No primary alcohol or phosphocholine group found"

    # Look for amino group patterns (free NH2, NH3+, or N-acylated)
    amino_patterns = [
        Chem.MolFromSmarts("[NH2,NH3+][CH]"), # Free amino
        Chem.MolFromSmarts("[NH]C(=O)"), # N-acylated
        Chem.MolFromSmarts("[NH+](C)C") # N,N-dimethyl
    ]
    
    has_amino = any(mol.HasSubstructMatch(pattern) for pattern in amino_patterns if pattern is not None)
    if not has_amino:
        return False, "No characteristic amino group found"

    # Look for characteristic 2-amino-1,3-diol or similar motif
    sphingoid_core_patterns = [
        Chem.MolFromSmarts("[CH2][CH]([NH2,NH3+,NH])[CH](O)"), # Basic sphingoid core
        Chem.MolFromSmarts("[CH2][CH]([NH2,NH3+,NH])[CH](=O)"), # 3-dehydro variant
        Chem.MolFromSmarts("[CH2][CH]([NH2,NH3+,NH])[CH](O)[CH](O)") # Phytosphingosine-like
    ]
    
    has_core = any(mol.HasSubstructMatch(pattern) for pattern in sphingoid_core_patterns if pattern is not None)
    if not has_core:
        return False, "No characteristic sphingoid core structure found"

    # Look for long carbon chain
    chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")  # At least 10 carbons in a chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain found"

    # Additional checks for common modifications
    has_double_bond = mol.HasSubstructMatch(Chem.MolFromSmarts("C=C"))
    has_glucose = mol.HasSubstructMatch(Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1"))
    has_acyl = mol.HasSubstructMatch(Chem.MolFromSmarts("NC(=O)"))
    
    # Build classification reason
    reason = "Contains sphingoid core structure with "
    features = []
    if has_double_bond:
        features.append("unsaturation")
    if has_glucose:
        features.append("glucosyl modification")
    if has_acyl:
        features.append("N-acylation")
    if len(features) > 0:
        reason += ", ".join(features)
    else:
        reason += "required functional groups"

    return True, reason