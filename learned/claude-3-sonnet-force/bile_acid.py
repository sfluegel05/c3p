"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_bile_acid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are characterized by a steroid backbone, presence of a carboxyl group,
    and at least one hydroxy group. Additional features like conjugation with glycine/taurine
    and additional rings or substituents are also considered.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@@]3([C@H](C[C@@H]4C[C@@H](C[C@]3([C@H]2[C@@]1([H])C)C)C)C)C"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"
    
    # Check for at least one hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("O")
    if len(mol.GetSubstructMatches(hydroxy_pattern)) < 1:
        return False, "No hydroxy group found"
    
    # Check for additional features (conjugation, rings, substituents)
    features = []
    
    # Conjugation with glycine or taurine
    conj_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if mol.HasSubstructMatch(conj_pattern):
        features.append("Conjugated with glycine/taurine")
    
    # Additional rings or substituents
    additional_pattern = Chem.MolFromSmarts("[r5,r6,r7]")
    if mol.HasSubstructMatch(additional_pattern):
        features.append("Additional rings or substituents present")
    
    # Scoring based on features
    score = len(features)
    if score >= 2:
        return True, f"Bile acid with following features: {', '.join(features)}"
    elif score == 1:
        return True, f"Bile acid with the following feature: {features[0]}"
    else:
        return True, "Bile acid with no additional features"