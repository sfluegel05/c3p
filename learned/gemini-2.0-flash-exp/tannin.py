"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are astringent polyphenolic vegetable principles, chiefly complex glucosides of catechol and pyrogallol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for multiple aromatic rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_matches) < 2:
         return False, "Too few aromatic rings for tannin"

    # 2. Check for hydroxyl groups (at least 4-5 is reasonable)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 4:
         return False, "Too few hydroxyl groups for tannin"
    
    # 3. Check for catechol or pyrogallol units (1,2 or 1,2,3 substituted benzene)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern) and not mol.HasSubstructMatch(pyrogallol_pattern):
        return False, "No catechol or pyrogallol units found."

    # 4. Check for glycosidic bonds (C-O-C with an adjacent 6-membered ring containing O - a simplified sugar check)
    glycosidic_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[C]1[CH2][CH2][CH][CH][CH][O]1")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    
    # 5. Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
       if len(glycosidic_matches) == 0: #If no sugar and mw is too low
           return False, "Molecular weight too low for tannin"
    if mol_wt > 2500:
        return False, "Molecular weight too high for a typical tannin"

    #If we reach this point, likely a tannin
    return True, "Molecule contains multiple aromatic rings, hydroxyl groups, and catechol/pyrogallol units, and has appropriate MW"