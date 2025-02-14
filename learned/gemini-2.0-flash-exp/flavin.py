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

    # Define the core dimethylisoalloxazine structure using a more specific SMARTS pattern
    core_pattern = Chem.MolFromSmarts("c1cc2nc3c(nc(=O)[nH]c3=O)nc2cc1")
    if not mol.HasSubstructMatch(core_pattern):
         return False, "Core dimethylisoalloxazine structure not found"

    # Get the nitrogen atom (N10) - this is the nitrogen directly bonded to the core.
    n10_match = mol.GetSubstructMatches(Chem.MolFromSmarts("n1c2cc(c)c(c)cc2nc2c1nc(=O)[nH]c2=O"))
    if not n10_match:
        return False, "N10 not found"
    
    # Check for a substituent at N10 (outside the core)
    for match in n10_match:
        n10_atom_index = -1
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "n":
                 n10_atom_index = atom_index
                 break

        if n10_atom_index == -1:
            return False, "N10 atom not found within substructure match."
        
        n10_atom = mol.GetAtomWithIdx(n10_atom_index)
        
        
        # Ensure that the degree is higher than 2, meaning it has a substituent outside the core
        if n10_atom.GetTotalDegree() > 2:
           return True, "Contains dimethylisoalloxazine core with a substituent on the 10 position"

    return False, "No substituent found at N10"