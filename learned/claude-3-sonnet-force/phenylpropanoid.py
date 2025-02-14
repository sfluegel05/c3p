"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26266 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenylpropane backbone (benzene ring with alkyl chain of 1-6 carbons)
    backbone_pattern = Chem.MolFromSmarts("[c]1ccccc1CCC[C,#0]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No phenylpropane backbone found"
    
    # Look for common substituents (hydroxy, methoxy, ester, prenyl, etc.) attached to the backbone
    substituents = ['O', 'OC', 'C(=O)O', 'C(=O)OC', 'c1ccccc1', 'CC=C(C)C']
    substituent_pattern = Chem.MolFromSmarts("|".join([f"[{sub}]~[#6]~[{backbone_pattern.GetSmarts()}]" for sub in substituents]))
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "No common phenylpropanoid substituents found attached to the backbone"
    
    # Look for additional structural features common to phenylpropanoids
    additional_patterns = [
        Chem.MolFromSmarts("[c]1ccc(cc1)C=CC(=O)c2ccccc2"),  # chalcones
        Chem.MolFromSmarts("[c]1ccc(cc1)C=Cc2ccccc2"),  # stilbenes
        Chem.MolFromSmarts("c1c(O)cc2c(c1)OC(C=3C(=O)OC(c3=O)c4ccccc4)CC2"),  # flavonoid scaffold
        Chem.MolFromSmarts("C1=C(C(=O)OC1)c2ccccc2"),  # coumarin scaffold
        Chem.MolFromSmarts("c1ccc(C[C@H](C)Cc2ccc(OC)cc2)cc1OC"),  # lignan scaffold
        Chem.MolFromSmarts("[O-][C+](=O)[C]~[c]"),  # ester or carboxylic acid
        Chem.MolFromSmarts("[O]C=CC=CC=C"),  # prenyl or related group
    ]
    
    for pattern in additional_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains phenylpropane backbone and additional structural features common to phenylpropanoids"
    
    return False, "Meets the basic requirements of a phenylpropanoid but lacks additional structural features"