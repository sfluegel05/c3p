"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
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
    
    score = 0 # Start with score of zero
    reason = ""

    # 1. Check for catechin core (more specific)
    catechin_core_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(C[C@H]2[C@H](O)[CH](O)c3c(O)cc(O)cc32)c1")
    catechin_core_matches = mol.GetSubstructMatches(catechin_core_pattern)
    score += len(catechin_core_matches) * 4  # High score for catechin core

    # 2. Check for inter-catechin linkages
    # C-C linkages between catechin core
    inter_catechin_cc_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(C[C@H]2[C@H](O)[CH](O)c3c(O)cc(O)cc32)c1-[CX4]-c4cc(O)c(O)c(C[C@H]5[C@H](O)[CH](O)c6c(O)cc(O)cc65)c4")
    inter_catechin_cc_matches = mol.GetSubstructMatches(inter_catechin_cc_pattern)
    score += len(inter_catechin_cc_matches) * 5
   
    # Ether linkages between catechin cores
    inter_catechin_ether_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(C[C@H]2[C@H](O)[CH](O)c3c(O)cc(O)cc32)c1-[OX2]-c4cc(O)c(O)c(C[C@H]5[C@H](O)[CH](O)c6c(O)cc(O)cc65)c4")
    inter_catechin_ether_matches = mol.GetSubstructMatches(inter_catechin_ether_pattern)
    score += len(inter_catechin_ether_matches) * 5

    # 3. Check for galloyl linkages to catechin
    galloyl_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(C(=O)O)c(O)c1-[OX2]-c2cc(O)c(O)c(C[C@H]3[C@H](O)[CH](O)c4c(O)cc(O)cc43)c2")
    galloyl_matches = mol.GetSubstructMatches(galloyl_pattern)
    score += len(galloyl_matches) * 4

    # 4. Check for gallic acid component
    gallic_acid_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(C(=O)O)c(O)c1")
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pattern)
    score += len(gallic_matches)


    # 5. Check for glycosidic bonds (more specific)
    glycosidic_pattern_6 = Chem.MolFromSmarts("[CX4]-[OX2]-[C]1[CH2][CH2][CH][CH][CH][O]1") #6 membered ring
    glycosidic_pattern_5 = Chem.MolFromSmarts("[CX4]-[OX2]-[C]1[CH2][CH2][CH][CH][O]1") #5 membered ring
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern_6) + mol.GetSubstructMatches(glycosidic_pattern_5)
    score += len(glycosidic_matches)

    # 6. Penalize flavonoid core structure
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(O)c(c2[o+]c(c(c(O)c2)c3cc(O)cc(O)c3)c1") # Basic Flavonoid pattern
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    if len(flavonoid_matches) > 0:
      score -= len(flavonoid_matches) * 3 # Penalize flavonoids.

    # 7. Molecular weight and size
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 500 <= mol_wt <= 3000:
      score += 3
    elif mol_wt > 3000: #Penalise too large
        score -= 2
    
    num_atoms = mol.GetNumAtoms()
    if num_atoms >= 30:
      score += 1
    
    # Final classification based on score.
    if score >= 12:
        return True, "Molecule contains multiple linked aromatic rings, hydroxyl groups, and catechol/pyrogallol units, with catechin core, galloyl groups, glycosidic bonds and appropriate MW"
    elif score >= 8:
         return False, "Molecule has some tannin-like features, but not enough"
    else:
        return False, "Molecule does not resemble a tannin"