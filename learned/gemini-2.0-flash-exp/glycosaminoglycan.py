"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are polysaccharides containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for presence of aminomonosaccharide rings using SMARTS patterns that include rings and nitrogens
    amino_sugar_pattern1 = Chem.MolFromSmarts("[C1][O][C]([C])[C][C][C]1[NX3]") # 6-membered ring with N
    amino_sugar_pattern2 = Chem.MolFromSmarts("[C1][O][C]([C])[C][C]1[NX3]")  # 5-membered ring with N
    amino_sugar_pattern3 = Chem.MolFromSmarts("[C1][O][C]([C])[C][C][C]1[NX3][CX3](=[OX1])") # 6-membered ring with N-acetyl
    amino_sugar_pattern4 = Chem.MolFromSmarts("[C1][O][C]([C])[C][C]1[NX3][CX3](=[OX1])")  # 5-membered ring with N-acetyl

    amino_sugar_matches1 = mol.GetSubstructMatches(amino_sugar_pattern1)
    amino_sugar_matches2 = mol.GetSubstructMatches(amino_sugar_pattern2)
    amino_sugar_matches3 = mol.GetSubstructMatches(amino_sugar_pattern3)
    amino_sugar_matches4 = mol.GetSubstructMatches(amino_sugar_pattern4)

    if not (amino_sugar_matches1 or amino_sugar_matches2 or amino_sugar_matches3 or amino_sugar_matches4):
        return False, "No aminomonosaccharide unit detected."

    # 2. Check for glycosidic linkages between two ring carbons.
    glycosidic_link_pattern = Chem.MolFromSmarts("[CX4;R]-[OX2]-[CX4;R]") # Carbon of ring linked to oxygen linked to another carbon ring

    glycosidic_matches = mol.GetSubstructMatches(glycosidic_link_pattern)
    
    if len(glycosidic_matches) == 0:
         return False, "No glycosidic linkage between rings detected."


    # 3. Minimal molecular weight check - GAGs are polysaccharides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a polysaccharide"


    return True, "Contains aminomonosaccharide units linked by glycosidic bonds."