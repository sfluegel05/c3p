from rdkit import Chem
from rdkit.Chem import AllChem

def is_triazolobenzodiazepine(smiles: str):
    """
    Determines if a molecule is a triazolobenzodiazepine (benzodiazepine ortho-fused with triazole).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a triazolobenzodiazepine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for presence of required rings
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() >= 3:
        return False, "Molecule does not contain at least 3 rings"

    # SMARTS pattern for triazolobenzodiazepine core structure
    # Matches benzodiazepine fused with triazole
    pattern = Chem.MolFromSmarts('c1ccc2c(c1)C(=NCn1ncnc1-2)c1ccccc1')
    
    if not mol.HasSubstructMatch(pattern):
        # Try alternative pattern with slightly different arrangement
        pattern2 = Chem.MolFromSmarts('c1ccc2c(c1)-n1c(C(=NC2)c2ccccc2)nnc1')
        if not mol.HasSubstructMatch(pattern2):
            # Try third pattern variation
            pattern3 = Chem.MolFromSmarts('c1ccc2c(c1)-n1c(C(=NC2)c2ccccc2Cl)nnc1')
            if not mol.HasSubstructMatch(pattern3):
                return False, "Does not contain triazolobenzodiazepine core structure"

    # Additional check for triazole ring
    triazole_pattern = Chem.MolFromSmarts('[#7]1[#7]=[#7][#6][#7]1')
    if not mol.HasSubstructMatch(triazole_pattern):
        return False, "Does not contain triazole ring"

    # Check for benzodiazepine pattern
    benzo_pattern = Chem.MolFromSmarts('c1ccc2c(c1)C(=NCn1ncnc1-2)c1ccccc1')
    if not mol.HasSubstructMatch(benzo_pattern):
        return False, "Does not contain benzodiazepine structure"

    return True, "Contains triazolobenzodiazepine core structure with proper ring fusion"
# Pr=None
# Recall=0.0