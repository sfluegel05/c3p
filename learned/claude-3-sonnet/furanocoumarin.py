"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:47835 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic requirements - must have oxygen atoms and multiple rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    o_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if o_count < 2:
        return False, "Insufficient number of oxygen atoms"

    # More comprehensive coumarin patterns
    coumarin_patterns = [
        # Basic coumarin core with variations in bond types and substitutions
        Chem.MolFromSmarts("O=C1Oc2ccccc2C=C1"),  # Basic coumarin
        Chem.MolFromSmarts("O=C1Oc2ccccc2CC1"),   # Dihydrocoumarin
        Chem.MolFromSmarts("O=C1Oc2cc[c,C]cc2C=C1"),  # Allow for substitutions
        Chem.MolFromSmarts("O=C1Oc2c([c,C])c([c,C])c([c,C])c2C([c,C])=C1"), # More general
        # Extended patterns for different fusion types
        Chem.MolFromSmarts("O=C1Oc2cc3c(cc2C=C1)"),  # For linear fusion
        Chem.MolFromSmarts("O=C1Oc2c3c(cc2C=C1)"),   # For angular fusion
    ]
    
    has_coumarin = any(mol.HasSubstructMatch(pat) for pat in coumarin_patterns if pat is not None)
    if not has_coumarin:
        return False, "No coumarin core found"

    # Comprehensive furan patterns
    furan_patterns = [
        Chem.MolFromSmarts("c1cocc1"),            # Basic furan
        Chem.MolFromSmarts("C1COC=C1"),           # Dihydrofuran
        Chem.MolFromSmarts("[c,C]1[c,C]o[c,C][c,C]1"),  # General furan
        Chem.MolFromSmarts("C1OC=CC1"),           # Alternative furan
        Chem.MolFromSmarts("C1=COC=C1")           # Another furan variant
    ]
    
    has_furan = any(mol.HasSubstructMatch(pat) for pat in furan_patterns if pat is not None)
    if not has_furan:
        return False, "No furan ring found"

    # Comprehensive fusion patterns for different types of furanocoumarins
    fusion_patterns = [
        # Linear furanocoumarins (furo[3,2-g] and furo[2,3-g])
        Chem.MolFromSmarts("O=C1Oc2cc3occc3cc2C=C1"),
        Chem.MolFromSmarts("O=C1Oc2cc3c(cc2C=C1)occ3"),
        # Angular furanocoumarins (furo[3,2-h] and furo[2,3-h])
        Chem.MolFromSmarts("O=C1Oc2c3occc3ccc2C=C1"),
        Chem.MolFromSmarts("O=C1Oc2c3c(occ3)ccc2C=C1"),
        # More general fusion patterns
        Chem.MolFromSmarts("O=C1Oc2c([c,C])c3oc[c,C]c3[c,C]c2C=C1"),
        Chem.MolFromSmarts("O=C1Oc2[c,C]c3oc[c,C]c3[c,C]c2C=C1"),
        # Patterns for dihydrofuranocoumarin variants
        Chem.MolFromSmarts("O=C1Oc2cc3C4COC=C4cc3cc2C=C1"),
        Chem.MolFromSmarts("O=C1Oc2c3C4COC=C4ccc3cc2C=C1")
    ]

    has_proper_fusion = any(mol.HasSubstructMatch(pat) for pat in fusion_patterns if pat is not None)
    if not has_proper_fusion:
        return False, "No proper fusion between furan and coumarin found"

    # Determine fusion type
    fusion_type = []
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in [
        "O=C1Oc2cc3occc3cc2C=C1",
        "O=C1Oc2cc3c(cc2C=C1)occ3"
    ]):
        fusion_type.append("linear")
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in [
        "O=C1Oc2c3occc3ccc2C=C1",
        "O=C1Oc2c3c(occ3)ccc2C=C1"
    ]):
        fusion_type.append("angular")

    fusion_desc = f"({', '.join(fusion_type)} fusion)" if fusion_type else "(fusion type undetermined)"
    return True, f"Furanocoumarin with {fusion_desc}"