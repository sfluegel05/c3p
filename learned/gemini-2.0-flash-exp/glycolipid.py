"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a lipid (fatty acids, diacylglycerol, or ceramide) linked to a carbohydrate (mono-, di-, or trisaccharide) through a glycosidic bond.
    Some glycolipids don't have glycerol or sphingosine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # 1. Basic Checks
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Glycosidic Linkage
    glycosidic_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic linkage (C-O-C) found"

    # 3. & 5. Sugar Rings
    sugar_ring_pattern = Chem.MolFromSmarts("[O][C][C][C][C][C]")
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    sugar_ring_pattern2 = Chem.MolFromSmarts("[O][C][C][C][C]")
    sugar_matches2 = mol.GetSubstructMatches(sugar_ring_pattern2)
    if not sugar_matches and not sugar_matches2:
         return False, "No sugar ring found"
    
    # 4. Lipid Component
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
         return False, "No fatty acid chains detected"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
         return False, "Chains too short to be fatty acids"

    # Check for diacylglycerol component 
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    
    # Check for ceramide component 
    ceramide_pattern = Chem.MolFromSmarts("[CX4][NX3][CX4]~[CX4]")
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)

    # If neither of these patterns is present and only one fatty acid chain was identified
    if not has_glycerol and not has_ceramide and len(fatty_acid_matches) == 1:
       pass
    #If none of the patterns are present return false
    elif not has_glycerol and not has_ceramide and len(fatty_acid_matches) > 1:
        return False, "No glycerol or ceramide found for more than one fatty acid"

    return True, "Glycosidic linkage, carbohydrate and lipid chains are present."