"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B Vitamins
Any member of the group of eight water-soluble vitamins originally thought to be a single compound (vitamin B) that play important roles in cell metabolism.
The group comprises of vitamin B1, B2, B3, B5, B6, B7, B9, and B12 (Around 20 other compounds were once thought to be B vitamins but are no longer classified as such).
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_B_vitamin(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Thiamine (B1)
    thiamine_patterns = ['*c1c([nH+]cs1)*', '*c1c([nH+]ccn1)*']
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in thiamine_patterns):
        return True, "Thiamine (B1) substructure found"

    # Riboflavin (B2)
    riboflavin_pattern = '*c1nc2c(=O)[nH]c(=O)[nH]c2c(*)cc1*'
    if mol.HasSubstructMatch(Chem.MolFromSmarts(riboflavin_pattern)):
        return True, "Riboflavin (B2) substructure found"

    # Niacin (B3)
    niacin_patterns = ['*c1ccncc1*', '*c1ccn(O)cc1*']
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in niacin_patterns):
        return True, "Niacin (B3) substructure found"

    # Pantothenic acid (B5)
    pantothenic_pattern = '*CC(C)(CO)C(O)C(=O)NCCC(=O)O*'
    if mol.HasSubstructMatch(Chem.MolFromSmarts(pantothenic_pattern)):
        return True, "Pantothenic acid (B5) substructure found"

    # Pyridoxine (B6)
    pyridoxine_patterns = ['*c1c(CO)ncc(CO)c1O*', '*c1c(CO)ncc(COP(O)(O)=O)c1O*']
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in pyridoxine_patterns):
        return True, "Pyridoxine (B6) substructure found"

    # Biotin (B7)
    biotin_pattern = '*C12CSC(CCCCC(O)=O)C1NC(=O)N2*'
    if mol.HasSubstructMatch(Chem.MolFromSmarts(biotin_pattern)):
        return True, "Biotin (B7) substructure found"

    # Folate (B9)
    folate_patterns = ['*c1cnc2nc(N)[nH]c(=O)c12*', '*c1cnc2nc(N)[nH]c(N)c12*']
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in folate_patterns):
        return True, "Folate (B9) substructure found"

    # Cobalamin (B12)
    cobalamin_pattern = '*[Co-3]1234*'
    if mol.HasSubstructMatch(Chem.MolFromSmarts(cobalamin_pattern)):
        return True, "Cobalamin (B12) substructure found"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rings = rdMolDescriptors.CalcNumRings(mol)

    if mol_wt < 100 or n_rings < 2:
        return False, "Molecular weight or ring count too low for B vitamin"

    return False, "No B vitamin substructures found"