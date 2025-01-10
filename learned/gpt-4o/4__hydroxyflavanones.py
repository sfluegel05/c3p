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

    # Define the core flavanone structure with a carbonyl group and a three-carbon bridge
    flavanone_core_smarts = "c1ccccc1C2C(=O)C3C=C(c4ccccc4)O3"
    flavanone_core_mol = Chem.MolFromSmarts(flayanone_core_smarts)
    if not mol.HasSubstructMatch(flayanone_core_mol):
        return False, "Missing core flavanone structure"
    
    # Pattern to specifically find a hydroxyl group at the 4' position on the B-ring
    hydroxy_b_ring_smarts = "c1cc(ccc1O)c2cc(c(=O)C3OC(C=C3)c2)c4ccccc4"  # Account for the correct hydroxyl position
    b_ring_hydroxy_mol = Chem.MolFromSmarts(hydroxy_b_ring_smarts)
    if mol.HasSubstructMatch(b_ring_hydroxy_mol):
        return True, "Molecule is a 4'-hydroxyflavanone"

    return False, "No 4'-hydroxy group found on flavanone B ring"