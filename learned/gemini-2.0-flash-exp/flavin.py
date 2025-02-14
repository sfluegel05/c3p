"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of dimethylisoalloxazine with a substituent on the 10 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core dimethylisoalloxazine structure using a SMARTS pattern
    core_pattern = Chem.MolFromSmarts("Cc1cc2Nc3c([nH]c(=O)[nH]c3=O)Nc2cc1C")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core dimethylisoalloxazine structure not found"

    # SMARTS pattern to specifically match N10 and its substituent
    # Map the N10 atom (atom 1), and any atom attached to N10 (atom 2)
    n10_pattern = Chem.MolFromSmarts("[#7:1]~[*:2]")

    # Get all matches for the N10 substructure
    n10_matches = mol.GetSubstructMatches(n10_pattern)

    # Check if at least one match exists that includes the correct N10 (index 0)
    found_n10_substituent = False
    for match in n10_matches:
        # The N10 in the core is the same N as in the pattern
        if match:
            core_mol = Chem.MolFromSmarts("c1cc2c(cc1C)N(C)C3=NC(=O)NC(=O)C3N2")
            core_match = mol.GetSubstructMatch(core_mol)
            if core_match:
                 # Get the correct N10 atom in the mol
                 core_n10_index = core_match[5]
                 
                 if match[0] == core_n10_index:
                     found_n10_substituent = True
                     break;
    
    if found_n10_substituent:
         return True, "Contains dimethylisoalloxazine core with a substituent on the 10 position"
    else:
        return False, "No substituent found at N10"