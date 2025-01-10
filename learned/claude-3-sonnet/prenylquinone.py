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
    
    # More comprehensive quinone patterns
    quinone_patterns = [
        "C1(=O)C=CC(=O)C=C1",  # para-quinone
        "C1(=O)C(=O)C=CC=C1",   # ortho-quinone
        "C1(=O)C(=O)C(*)=C(*)C(*)=C1*",  # substituted quinone
        "C1(=O)C(*)=C(*)C(=O)C(*)=C1*",  # substituted para-quinone
        "O=C1C(=C)C(=O)C=CC1",  # alternative quinone
        "O=C1C=CC(=O)C(*)=C1*"  # substituted quinone variant
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
        "CC(C)=CCCC(C)=C",  # diprenyl
        "CC(C)=CCCC(=C)C",  # alternative prenyl
        "C/C=C(/C)CC/C=C(/C)C",  # trans-prenyl chain
        "CC(C)=CCCC(C)CC=C"  # longer prenyl variant
    ]
    
    has_prenyl = False
    for pattern in prenyl_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_prenyl = True
            break
            
    if not has_prenyl:
        return False, "No prenyl/polyprenyl chain found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for a prenylquinone"
        
    if o_count < 2:
        return False, "Too few oxygens for a quinone structure"

    # Improved subclass patterns
    menaquinone_pattern = Chem.MolFromSmarts("[#6]1=CC=CC2=C1C(=O)C([CH3,CH2])=C([CH2,CH])C2=O")
    ubiquinone_pattern = Chem.MolFromSmarts("COC1=C(OC)C(=O)C([CH2,CH])=C([CH3,CH2])C1=O")
    plastoquinone_pattern = Chem.MolFromSmarts("CC1=C(C)C(=O)C=C(CC=C)C1=O")
    
    # Calculate rotatable bonds and ring count
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Check for characteristic substitution patterns
    methoxy_pattern = Chem.MolFromSmarts("COC")
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    
    # Identify specific subclasses
    features = []
    if mol.HasSubstructMatch(menaquinone_pattern):
        features.append("menaquinone-like")
    if mol.HasSubstructMatch(ubiquinone_pattern):
        features.append("ubiquinone-like")
    if mol.HasSubstructMatch(plastoquinone_pattern):
        features.append("plastoquinone-like")
    if has_methoxy:
        features.append("methoxy-substituted")
    if n_rotatable >= 7:
        features.append("long prenyl chain")
        
    # Stricter classification criteria
    if len(features) >= 1 and has_prenyl and has_quinone:
        return True, f"Prenylquinone with features: {', '.join(features)}"
    elif has_prenyl and has_quinone and n_rotatable >= 5:
        return True, "Basic prenylquinone structure with prenyl chain"
    else:
        return False, "Missing characteristic prenylquinone features"