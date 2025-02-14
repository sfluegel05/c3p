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

    # Define the core dimethylisoalloxazine structure using a more specific SMARTS pattern with atom mapping
    # Match the core structure and *map* the N10 nitrogen (atom 1) and the connected atom (atom 2)
    core_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]2[#7][#6]3[#7](-[#6](=[#8])[#7]3=[#8])[#7]2[#6][#6]1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core dimethylisoalloxazine structure not found"

    # SMARTS pattern to specifically match N10 and its substituent
    # The key is mapping the N10 atom (atom 1), and any atom attached to N10 (atom 2) which represents the substituent
    n10_pattern = Chem.MolFromSmarts("[#7:1]~[*:2]")
    
    n10_matches = mol.GetSubstructMatches(core_pattern)
    
    found_n10_substituent = False

    for core_match in n10_matches:
        #Get the nitrogen corresponding to the core from the first pattern
         core_n10_index = -1
         for atom_index in core_match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "n" and atom.GetTotalDegree() == 3:
               core_n10_index = atom_index
               break
        if core_n10_index == -1:
            return False, "Core N10 not found"


        n10_substituent_matches = mol.GetSubstructMatches(n10_pattern)
        
        for sub_match in n10_substituent_matches:
                if len(sub_match) == 2:
                    mapped_n10_index = sub_match[0]
                    if mapped_n10_index == core_n10_index:
                        found_n10_substituent = True
                        break
        if found_n10_substituent:
                break
        

    if found_n10_substituent:
         return True, "Contains dimethylisoalloxazine core with a substituent on the 10 position"
    else:
        return False, "No substituent found at N10"