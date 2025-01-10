"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined as a hexose with the 'D' configuration at position 5.

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

    # Check if it's a hexose (6 carbons, typically two types of rings)
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons != 6:
        return False, "Not a hexose, incorrect number of carbon atoms"
    
    # Check for typical D-hexose ring structures
    hexose_ring_patterns = [
        Chem.MolFromSmarts('OC1O[C@@H](O)C[C@@H](O)CO1'),  # pyranose form
        Chem.MolFromSmarts('OC1OC[C@@H](O)C[C@@H]1O')      # furanose form
    ]
    
    is_hexose_ring = any(mol.HasSubstructMatch(pattern) for pattern in hexose_ring_patterns)
    if not is_hexose_ring:
        return False, "No hexose ring structure found"
    
    # Check stereochemistry: D-configuration typically specific stereo at carbon 5, adjust for typical representation here
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    d_stereo_count = 0
    for center in stereo_centers:
        atom_idx, chirality = center
        atom = mol.GetAtomWithIdx(atom_idx)
        # Assume position of carbon-5 in stereo centers as an application case, might need further specific mapping
        if chirality == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:  # or some further specific check
            d_stereo_count += 1
    
    if d_stereo_count == 0:
        return False, "No D-configuration stereochemistry at position 5 detected"
        
    return True, "Molecule identified as D-hexose with correct stereochemistry"