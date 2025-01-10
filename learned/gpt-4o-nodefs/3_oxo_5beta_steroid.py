"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid should have a steroid backbone with a ketone group at the 3-position
    and specific stereochemistry denoted as '5beta'.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more comprehensive steroid backbone pattern
    steroid_smarts = "[#6]1:[#6][#6]2[#6][#6]3[#6][#6]4[#6,#8]([#6]3[#6]2[#6]4)[#6]1"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"
    
    # Define a pattern for the 3-oxo group specifically on the steroid framework
    oxo_pattern = Chem.MolFromSmarts("[$([#6;R]=[O;R])]")  # Ring-embedded C=O
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found on steroid framework"
    
    # Evaluate stereochemistry for 5beta orientation
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    beta_orientation = False
    for atom_idx, chirality in chiral_centers:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic() or atom.GetTotalNumHs() == 0:
            continue
        if "S" in chirality:  # Check if the chirality suggests beta orientation
            beta_orientation = True
            break
    
    if not beta_orientation:
        return False, "5beta stereochemistry not confirmed"
    
    return True, "Molecule matches the 3-oxo-5beta-steroid characteristics"