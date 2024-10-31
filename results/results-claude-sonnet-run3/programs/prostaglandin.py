from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for cyclopentane ring which is characteristic of prostaglandins
    rings = mol.GetRingInfo()
    has_5_ring = False
    for ring in rings.AtomRings():
        if len(ring) == 5:
            has_5_ring = True
            break
            
    if not has_5_ring:
        return False, "No cyclopentane ring found"
        
    # Identify specific prostaglandin type based on characteristic features
    pga_pattern = Chem.MolFromSmarts('[*][C@@H]1[C@@H]([*])C=CC1=O')
    pgb_pattern = Chem.MolFromSmarts('[*]C1=C([*])CCC1=O')
    pgc_pattern = Chem.MolFromSmarts('[*][C@H]1C(=O)CC=C1[*]')
    pgd_pattern = Chem.MolFromSmarts('O[C@H]1CC(=O)[C@H]([*])[C@H]1[*]')
    pge_pattern = Chem.MolFromSmarts('O[C@@H]1CC(=O)[C@H]([*])[C@H]1[*]')
    pgf_pattern = Chem.MolFromSmarts('OC1C[C@@H](O)[C@H]([*])[C@H]1[*]')
    pgg_pattern = Chem.MolFromSmarts('[*][C@H]1[C@@H]2C[C@@H](OO2)[C@@H]1[*]')
    pgh_pattern = Chem.MolFromSmarts('[*][C@H]1[C@@H]2C[C@@H](OO2)[C@@H]1[*]')
    pgi_pattern = Chem.MolFromSmarts('[H][C@]12C[C@@H](O)[C@H]([*])[C@@]1([H])CC([*])O2')
    pgj_pattern = Chem.MolFromSmarts('[*][C@H]1C=CC(=O)[C@@H]1[*]')

    if mol.HasSubstructMatch(pga_pattern):
        return True, "Prostaglandin A"
    elif mol.HasSubstructMatch(pgb_pattern):
        return True, "Prostaglandin B"
    elif mol.HasSubstructMatch(pgc_pattern):
        return True, "Prostaglandin C"
    elif mol.HasSubstructMatch(pgd_pattern):
        return True, "Prostaglandin D"
    elif mol.HasSubstructMatch(pge_pattern):
        return True, "Prostaglandin E"
    elif mol.HasSubstructMatch(pgf_pattern):
        return True, "Prostaglandin F"
    elif mol.HasSubstructMatch(pgg_pattern):
        return True, "Prostaglandin G"
    elif mol.HasSubstructMatch(pgh_pattern):
        return True, "Prostaglandin H"
    elif mol.HasSubstructMatch(pgi_pattern):
        return True, "Prostaglandin I"
    elif mol.HasSubstructMatch(pgj_pattern):
        return True, "Prostaglandin J"
        
    return False, "Does not match prostaglandin structural patterns"
# Pr=1.0
# Recall=1.0