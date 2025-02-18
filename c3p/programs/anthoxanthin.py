"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments, with a benzopyran-4-one core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core benzopyran-4-one structure (flavone or flavonol) with carbonyl group
    # More robust SMARTS pattern that captures flavones, flavonols, and isoflavones
    core_pattern = Chem.MolFromSmarts("[c]1[c]([c][c][c]2)[o][c]([c](=[O])[c]3[c][c][c][c][c]3)[c]2[c]1 | [c]1[c]([c]2[c]([c](=[O])[c]3[c][c][c][c][c]3)[o][c]2[c]1)c1ccccc1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No benzopyran-4-one core structure found"

    #2. Check for multiple rings.
    if rdMolDescriptors.CalcNumRings(mol) < 2:
        return False, "Molecule has less than 2 rings, not a flavonoid"

    # Check for a carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Molecule does not have a carbonyl group"

    # 3. Check for at least 2 hydroxyl groups.
    num_hydroxyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    num_methoxys = len(mol.GetSubstructMatches(Chem.MolFromSmarts("OC")))
    if num_hydroxyls < 2 :
        return False, "Less than 2 hydroxyl groups found"
    
    #4. Exclude molecules with glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][C]([C][O][C])([C][O][C])")
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Molecule has glycosidic bonds"

    #5 Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for an anthoxanthin"

    return True, f"Benzopyran-4-one core structure with {num_hydroxyls} hydroxyl and {num_methoxys} methoxy group(s) found"