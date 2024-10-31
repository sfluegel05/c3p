from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanin cation (a glycoside derivative of an anthocyanidin cation).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanin cation, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of positive charge
    has_positive_charge = False
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            has_positive_charge = True
            break
    
    if not has_positive_charge:
        return False, "No positive charge found - anthocyanins are cations"

    # Check for benzopyrylium core structure (characteristic of anthocyanins)
    # This SMARTS pattern represents the core structure with the O+ in the heterocyclic ring
    benzopyrylium_pattern = Chem.MolFromSmarts("c1c([O+]=c2cc(O)cc(O)c2)c(O)cc(O)c1")
    if not mol.HasSubstructMatch(benzopyrylium_pattern):
        benzopyrylium_pattern2 = Chem.MolFromSmarts("[O+]1=CC=C2C=C(O)C(O)=CC2=C1")
        if not mol.HasSubstructMatch(benzopyrylium_pattern2):
            return False, "No benzopyrylium core structure found"

    # Check for glycoside - look for characteristic sugar pattern
    # This looks for a pyranose ring connected through oxygen
    glycoside_patterns = [
        Chem.MolFromSmarts("OCC1OC(O)C(O)C(O)C1O"),  # Basic pyranose
        Chem.MolFromSmarts("OCC1OC(OC)C(O)C(O)C1O"),  # Modified pyranose
        Chem.MolFromSmarts("CC1OC(OC)C(O)C(O)C1O")    # Deoxy sugar
    ]
    
    has_glycoside = False
    for pattern in glycoside_patterns:
        if mol.HasSubstructMatch(pattern):
            has_glycoside = True
            break
            
    if not has_glycoside:
        return False, "No glycoside moiety found"

    # Check for characteristic substitution pattern
    # Anthocyanins typically have hydroxyl groups at specific positions
    hydroxyl_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic hydroxyl substitution pattern"

    # If all checks pass, this is likely an anthocyanin cation
    return True, "Molecule contains benzopyrylium core, glycoside moiety, and characteristic substitution pattern"
# Pr=None
# Recall=0.0