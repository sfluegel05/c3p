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
        "[#6]1(=O)[#6]=,:[#6][#6](=O)[#6]=,:[#6]1",  # basic quinone
        "O=C1C(=O)C=CC=C1",  # ortho-quinone
        "O=C1C=CC(=O)C=C1",  # para-quinone
        "O=C1C(O)=C(*)C(=O)C=C1*",  # hydroxyquinone
        "O=C1C(O)=C(*)C(=O)C(*)=C1*",  # substituted hydroxyquinone
        "O=C1C(*)=C(*)C(=O)C(*)=C1*",  # fully substituted quinone
        "O=C1[#6]=C([#6])[#6](=O)[#6]=C1",  # generic quinone core
        "O=C1C(=C)C(=O)C(*)C=C1",  # alternative quinone
    ]

    # More specific prenyl patterns
    prenyl_patterns = [
        "CC(C)=CCC[C@H]?(C)",  # basic prenyl with optional stereochem
        "C/C=C(/C)CC/C=C(/C)",  # extended prenyl
        "CC(C)=CCCC(=C)C",  # branched prenyl
        "C/C=C(/C)CC/C=C(/C)CC/C=C(/C)",  # long prenyl chain
        "[CH3][C](=C)[CH2][CH2][C](=C)[CH3]",  # isoprene units
    ]

    # Check for quinone core
    has_quinone = False
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_quinone = True
            break

    # Count prenyl units and check chain length
    prenyl_count = 0
    long_chain = False
    for pattern in prenyl_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        prenyl_count += len(matches)
        if len(matches) > 0 and rdMolDescriptors.CalcNumRotatableBonds(mol) >= 8:
            long_chain = True

    # Specific prenylquinone class patterns
    class_patterns = {
        "ubiquinone": "COC1=C(OC)C(=O)C(C/C=C(/C)C)=C(C)C1=O",
        "menaquinone": "CC1=C(C)C(=O)c2ccccc2C1=O",
        "plastoquinone": "CC1=C(C)C(=O)C=C(CC=C(C)C)C1=O",
        "tocopherolquinone": "CC1=C(O)C(=O)C(C)=C(C)C1=O"
    }

    # Check for specific class matches
    class_matches = []
    for class_name, pattern in class_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            class_matches.append(class_name)

    # Calculate descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Classification logic
    if not has_quinone:
        return False, "No quinone core structure found"

    if prenyl_count == 0 and not long_chain:
        return False, "No significant prenyl/polyprenyl chain found"

    # Strong evidence for prenylquinone
    if has_quinone and (long_chain or prenyl_count >= 2) and class_matches:
        return True, f"Matches known prenylquinone class: {', '.join(class_matches)}"

    # Check for characteristic features
    if has_quinone and prenyl_count >= 1:
        features = []
        if long_chain:
            features.append("long prenyl chain")
        if mol.HasSubstructMatch(Chem.MolFromSmarts("COC")):
            features.append("methoxy-substituted")
        if o_count >= 4:
            features.append("multiple oxygen functionalities")
        
        # Stricter criteria for classification
        if (long_chain or prenyl_count >= 2) and len(features) >= 2:
            if mol_wt >= 250 and c_count >= 15:
                return True, f"Prenylquinone with features: {', '.join(features)}"

    # Basic prenylquinone structure
    if has_quinone and long_chain and mol_wt >= 250 and n_rotatable >= 8:
        return True, "Basic prenylquinone structure with significant prenyl chain"

    return False, "Does not meet prenylquinone structural requirements"