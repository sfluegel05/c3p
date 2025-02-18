"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid (CHEBI:33525)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid has a carbohydrate linked via glycosidic bond to a lipid moiety,
    which can be a glycerol with acyl groups or a sphingoid base with fatty acid.
    Bacterial variants may have acylated carbohydrates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Improved carbohydrate detection: look for at least 3 contiguous hydroxyl-bearing carbons in a ring
    pyranose_pattern = Chem.MolFromSmarts("[C;R5][C;R5][C;R5][C;R5][C;R5][C;R5]([OH])")  # 6-membered with OH
    furanose_pattern = Chem.MolFromSmarts("[C;R5][C;R5][C;R5][C;R5][C;R5][OH]")  # 5-membered with OH
    if not (mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)):
        return False, "No carbohydrate moiety detected"
    
    # Detect glycosidic bond: oxygen connecting ring (carb) to non-ring atom (lipid part)
    glycosidic_pattern = Chem.MolFromSmarts("[C;R][OX2][!C;!R]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage between carb and lipid"
    
    # Check lipid part: either glycerol with acyl groups or sphingoid base with fatty acid
    # Glycerol case: 1,2-diacylglycerol
    glycerol = Chem.MolFromSmarts("[CH2X4][CHX4]([CH2X4])")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3]=[OX1]")
    if mol.HasSubstructMatch(glycerol):
        # Check for at least two ester groups attached to glycerol
        glycerol_esters = Chem.MolFromSmarts("[CH2X4][CHX4]([CH2X4])-([OX2][CX3]=[OX1])")
        if len(mol.GetSubstructMatches(glycerol_esters)) >= 2:
            return True, "1,2-Diacylglycerol with carbohydrate linkage"
    
    # Sphingoid base case: long chain with amine and hydroxyls
    sphingosine_pattern = Chem.MolFromSmarts("[CH2][CH]([NH2])[CH]([OH])[CH2]")
    if mol.HasSubstructMatch(sphingosine_pattern):
        # Check for amide-linked fatty acid
        amide_pattern = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
        if mol.HasSubstructMatches(amide_pattern):
            return True, "Sphingosine base with fatty acid and carbohydrate"
    
    # Bacterial case: acylated carbohydrate (at least one ester/amide on carb)
    carb_acyl = Chem.MolFromSmarts("[C;R][OX2][CX3]=[OX1]")  # Ester on carbohydrate ring
    if len(mol.GetSubstructMatches(carb_acyl)) >= 1:
        # Verify lipid chain length (at least 8 carbons in acyl group)
        for match in mol.GetSubstructMatches(carb_acyl):
            acyl_group = mol.GetAtomWithIdx(match[2]).GetNeighbors()[0]
            chain_length = rdMolDescriptors.CalcNumAdjacentHeavyAtoms(acyl_group)
            if chain_length >= 8:
                return True, "Acylated carbohydrate with lipid chain"
    
    return False, "Does not meet glycolipid structural criteria"