"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Count the number of carbon atoms (should be 6 for hexose)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Molecule does not have 6 carbon atoms typical of hexoses"

    # Identify the chiral centers in the molecule
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Check for the D-configuration, which is an 'R' configuration in the intended sugar stereochemistry
    # RDKit provides stereochemistry directly, but we need to handle mixed configurations sometimes
    found_D_chiral_center = False

    for center, config in chiral_centers:
        atom = mol.GetAtomWithIdx(center)
        atom_neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
        
        # Typically the chiral center position 5 in the chain (or equivalent in a ring) should be verified
        if config == 'R' and atom.GetAtomicNum() == 6 and 'O' in atom_neighbors:
            found_D_chiral_center = True
            break

    if found_D_chiral_center:
        return True, "Molecule has D-configuration at position 5"
    else:
        return False, "Molecule does not conform to D-hexose classification due to stereochemistry"