"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids typically contain around 20 carbons, exhibit complex ring structures,
    and a variety of functional groups.
    
    They often include polycyclic systems and can have varying functionalizations.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms - allowing for more flexibility
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15 or num_carbons > 40:
        return False, f"Contains {num_carbons} carbon atoms, outside typical diterpenoid range"

    # Check for multiple polycyclic systems
    num_rings = Chem.CalcNumRings(mol)
    if num_rings < 3:
        return False, "Diterpenoids typically have complex polycyclic ring systems"
    
    # Check for a broader set of functional groups
    functional_groups = [
        Chem.MolFromSmarts("[OX2H]"),  # Alcohol
        Chem.MolFromSmarts("C=O"),    # Ketone
        Chem.MolFromSmarts("[OX2]([#6])[#6]"),  # Ether
        Chem.MolFromSmarts("[NX1]#[CX2]"),   # Isocyano
        Chem.MolFromSmarts("COC(=O)"),  # Ester
        Chem.MolFromSmarts("C=C"),     # Olefin (Double bond)
        Chem.MolFromSmarts("OC=O"),    # Carboxylic acid
        Chem.MolFromSmarts("NC=O"),    # Amide
    ]
    
    has_functional_group = any(mol.HasSubstructMatch(pattern) for pattern in functional_groups)
    
    if not has_functional_group:
        return False, "Lacks common functional groups characterizing diterpenoids"
    
    return True, "Contains characteristics typical of diterpenoids, such as polycyclic ring structures and diverse functional groups"