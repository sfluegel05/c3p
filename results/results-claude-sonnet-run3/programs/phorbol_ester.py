from rdkit import Chem
from rdkit.Chem import AllChem

def is_phorbol_ester(smiles: str):
    """
    Determines if a molecule is a phorbol ester.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phorbol ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define SMARTS pattern for phorbol core structure
    # This pattern represents the characteristic tetracyclic ring system of phorbol with specific stereochemistry
    phorbol_core = Chem.MolFromSmarts('[C]1[C]2([H])[C]=C(CO)[C][C]3(O)[C](=O)[C]([CH3])=[C][C]3([H])[C]2(O)[CH]([CH3])[CH]([OH,O])[C]4[O,C][C]1[C]4(C)C')
    
    if not mol.HasSubstructMatch(phorbol_core):
        # Try alternative SMARTS pattern that might match different representations
        alt_core = Chem.MolFromSmarts('[C]12[C][C](CO)[C][C]3(O)[C](=O)[C]([CH3])=[C][C]3[C]1(O)[CH]([CH3])[CH]([OH,O])[C]1[O][C]2[C]1([CH3])[CH3]')
        if not mol.HasSubstructMatch(alt_core):
            return False, "Does not match phorbol core structure"

    # Check for presence of at least one ester group
    ester_pattern = Chem.MolFromSmarts('[#6]-C(=O)-O-[#6]')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester groups found"

    # Count number of ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_esters = len(ester_matches)

    # Additional structural checks
    # Check for required ketone at C-3 position
    ketone_pattern = Chem.MolFromSmarts('[#6]-C(=O)-[#6]')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing required ketone group"

    # Check for required hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 1:
        return False, "Insufficient number of hydroxyl groups"

    return True, f"Phorbol ester with {num_esters} ester group(s)"
# Pr=None
# Recall=0.0