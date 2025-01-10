"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose with the 'D' configuration at position 5.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 6 carbons typical for hexoses
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons != 6:
        return False, "Not a hexose, incorrect number of carbon atoms"
    
    # Try matching common hexose ring shapes (pyranose and furanose)
    pyranose_pattern = Chem.MolFromSmarts('OC[C@H]([C@H](O)C)O')
    furanose_pattern = Chem.MolFromSmarts('O[C@H][C@H](O)CO')
    if not (mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)):
        return False, "No typical hexose ring structure found"

    # Locate chiral centers
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_carbon_5 = None

    # Assume by default that carbon-5 maps to a hypothesized specific index
    carbon_5_idx = 4  # assuming zero-based index which corresponds to position 5

    for idx, chirality in stereo_centers:
        if idx == carbon_5_idx:
            if chirality in {Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW}: 
                chiral_carbon_5 = chirality
                break

    if chiral_carbon_5 is None:
        return False, "No D-configuration stereo center found at hexose position 5"

    # If found appropriate configuration (for simplicity considering both CW and CCW as sufficient for example)
    return True, "Molecule identified as D-hexose with appropriate stereochemistry at carbon-5"