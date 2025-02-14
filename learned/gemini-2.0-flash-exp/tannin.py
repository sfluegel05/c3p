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
    
    reason = ""
    
    # 1. Check for multiple aromatic rings (at least 3). This is a broader check
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1") # Basic aromatic ring
    aromatic_matches = mol.GetSubstructMatches(aromatic_pattern)
    if len(aromatic_matches) < 3:
        return False, "Fails: fewer than 3 aromatic rings."

    # 2. Check for multiple hydroxyl groups (at least 5)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 5:
        return False, f"Fails: fewer than 5 hydroxyl groups, got {len(hydroxyl_matches)}."

    # 3. Check for ester or ether linkages
    ester_ether_pattern = Chem.MolFromSmarts("[OX2]-[CX4]")
    ester_ether_matches = mol.GetSubstructMatches(ester_ether_pattern)
    if len(ester_ether_matches) < 3:
          return False, f"Fails: fewer than 3 ether/ester linkages, got {len(ester_ether_matches)}."
    
    # 4. Check for glycosidic bonds (C-O-C where one C is part of a sugar)
    glycosidic_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[C]1[CH2][CH2][CH](O)[CH][CH][O]1")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)


    # 5. Check for gallic acid and catechin units
    gallic_acid_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(C(=O)O)c(O)c1")
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pattern)
    catechin_core_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(C[C@H]2[C@H](O)[CH](O)c3c(O)cc(O)cc32)c1")
    catechin_core_matches = mol.GetSubstructMatches(catechin_core_pattern)

    if len(gallic_matches) < 1 and len(catechin_core_matches) <1 and len(glycosidic_matches) <1:
        return False, "Fails: does not contain catechin, gallic acid or glycosidic bonds"


    # 6. Molecular weight and atom counts (crude size check)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Fails: Molecular weight {mol_wt} is too low."
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 25:
         return False, f"Fails: Too few atoms, only {num_atoms}"

    # If all checks pass, then its a tannin
    return True, "Molecule contains multiple aromatic rings, hydroxyl groups, ether/ester linkages, and catechin/gallic acid/glycosidic substructures. Molecular weight and atom count are within range."