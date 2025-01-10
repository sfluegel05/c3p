"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined as a hexose with D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Molecule does not have 6 carbon atoms typical of hexoses"

    # Identify the chiral centers in the molecule
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Check for the D-configuration, which usually occurs at the 5th carbon in Fischer projection
    # Let's assume the SMILES provides the structure starting with the aldehyde end for open chain
    # Or the anomeric carbon for cyclic forms, we want the 5th C atom to have R configuration
    # Typically (R) stereochemistry implies 'D' configuration in sugars because of the convention
    # Here, verify chiral centers assuming a simple straight chain 6-carbon structure

    for center, config in chiral_centers:
        if config == 'R':
            carbon_idx = center + 1  # rdkit index is 0-based
            if carbon_idx == 5:
                # Checking for 5th carbon is R-type
                return True, "Molecule has D-configuration at position 5"
                
    return False, "Molecule does not conform to D-hexose classification due to stereochemistry"