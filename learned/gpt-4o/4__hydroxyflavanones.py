"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4_hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flavanone core with proper stereochemistry taking into account variable stereocenters
    flavanone_core_smarts = "[cH1]1[cH1][cH1][cH1][cH1][cH1]2[C@@H]([C@@H]([C@@H]2C1=O)[cH1]3[cH1][cH1][cH1][cH1][cH1]3)O" 
    flavanone_core_mol = Chem.MolFromSmarts(flayanone_core_smarts)
    if not mol.HasSubstructMatch(flayanone_core_mol):
        return False, "Missing core flavanone structure"
    
    # Correct pattern to specifically find hydroxyl attached at the 4' position
    hydroxy_b_ring_smarts = "O[cH]1[cH][cH][cH](O)[cH][cH]1"  # Apply general hydroxy to benzene using SMARTS
    b_ring_hydroxy_mol = Chem.MolFromSmarts(hydroxy_b_ring_smarts)
    if mol.HasSubstructMatch(b_ring_hydroxy_mol):
        return True, "Molecule is a 4'-hydroxyflavanone"

    return False, "No 4'-hydroxy group found on flavanone B ring"