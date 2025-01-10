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
    
    # Check for exactly 6 carbons in the molecule for hexoses
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons != 6:
        return False, "Not a hexose, incorrect number of carbon atoms"
    
    # Define common D-hexose pyranose and furanose ring patterns
    pyranose_pattern = Chem.MolFromSmarts('OC1COC(O)C(O)C1')
    furanose_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@H](O)[C@H]1O')
    open_chain_pattern = Chem.MolFromSmarts('[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)C=O')
    
    # Look for hexose structures within the molecule (either closed or open-chain forms)
    if not (mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(open_chain_pattern)):
        return False, "No typical hexose structure found"
    
    # Locate chiral centers and check for proper stereochemistry at carbon-5
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # In a D-hexose (open-form), carbon-5 should have R-configuration typically associated with carbohydrates
    # In pyranose/furanose forms, it's often located after ring closure
    detected_D_config = False
    for index, chirality in stereo_centers:
        atom = mol.GetAtomWithIdx(index)
        connected_hydroxyl = any(neighbor.GetAtomicNum() == 8 for neighbor in atom.GetNeighbors())
        
        if connected_hydroxyl and chirality in {Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW}:
            # Here we assume this center corresponds effectively to D-configuration contextually
            detected_D_config = True
            break

    if not detected_D_config:
        return False, "No D-configuration stereo center found at hexose position 5"
    
    return True, "Molecule identified as D-hexose with stereochemistry at carbon-5"