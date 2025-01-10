"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
A quinone substituted by a polyprenyl-derived side-chain.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded quinone patterns to catch more variants
    quinone_patterns = [
        "C1(=O)C=CC(=O)C=C1",  # para-quinone
        "C1(=O)C(=O)C=CC=C1",   # ortho-quinone
        "O=C1C(=C)C(=O)C=CC1",  # alternative quinone
        "O=C1C=CC(=O)C(*)=C1",  # substituted quinone
        "O=C1C(*)=C(*)C(=O)C(*)=C1*",  # heavily substituted quinone
        "O=C1C(O)=C(*)C(=O)C(*)=C1*",  # hydroxyquinone
        "[#6]1(=O)[#6]=,:[#6][#6](=O)[#6]=,:[#6]1",  # generic quinone
        "O=C1C(=C)C(=O)C(*)C(*)=C1*"  # complex quinone
    ]
    
    has_quinone = False
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_quinone = True
            break
            
    if not has_quinone:
        return False, "No quinone core structure found"

    # Improved prenyl/polyprenyl patterns
    prenyl_patterns = [
        "CC(C)=CCC",  # basic prenyl
        "C/C=C(/C)CC",  # trans prenyl
        "CC(C)=CCCC(C)=C",  # diprenyl
        "[CH3][C]=[C][CH2][CH2]",  # isoprene unit
        "C/C=C(/C)CC/C=C(/C)",  # extended prenyl
        "[CH3][C](=[CH2])[CH2][CH2]",  # terminal prenyl
        "CC(C)=CCCC(=C)C",  # branched prenyl
    ]
    
    prenyl_count = 0
    for pattern in prenyl_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        prenyl_count += len(matches)

    if prenyl_count == 0:
        return False, "No prenyl/polyprenyl chain found"

    # Calculate molecular descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Specific subclass patterns
    menaquinone_pattern = Chem.MolFromSmarts("[#6]1=CC=CC2=C1C(=O)C=C(CC=C)C2=O")
    ubiquinone_pattern = Chem.MolFromSmarts("COC1=C(OC)C(=O)C=C(CC=C)C1=O")
    plastoquinone_pattern = Chem.MolFromSmarts("CC1=C(C)C(=O)C=C(CC=C)C1=O")

    features = []
    if mol.HasSubstructMatch(menaquinone_pattern):
        features.append("menaquinone-like")
    if mol.HasSubstructMatch(ubiquinone_pattern):
        features.append("ubiquinone-like")
    if mol.HasSubstructMatch(plastoquinone_pattern):
        features.append("plastoquinone-like")
    
    # Check for characteristic substitution
    if mol.HasSubstructMatch(Chem.MolFromSmarts("COC")):
        features.append("methoxy-substituted")
    if n_rotatable >= 7:
        features.append("long prenyl chain")

    # Classification criteria
    is_prenylquinone = False
    reason = "Missing characteristic prenylquinone features"

    # Main classification logic
    if has_quinone and prenyl_count >= 1:
        if len(features) >= 1:
            is_prenylquinone = True
            reason = f"Prenylquinone with features: {', '.join(features)}"
        elif n_rotatable >= 5 and c_count >= 15 and o_count >= 2:
            # Check for minimum size and composition
            if mol_wt >= 200:  # Minimum weight for a basic prenylquinone
                is_prenylquinone = True
                reason = "Basic prenylquinone structure with prenyl chain"

    return is_prenylquinone, reason