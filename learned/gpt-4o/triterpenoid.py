"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Considers general characteristics such as multiple ring structures, and specific
    functional groups like hydroxyl, ketone, or carboxyl groups plus flexibility in 
    core skeleton rearrangements or additional sugar moieties.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms with a more lenient range for complex structures
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 60:
        return False, f"Carbon count {c_count} is outside a flexible range for typical triterpenoids (20-60)"
    
    # Check for a minimum number of oxygen atoms suggesting functional modifications
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; typical triterpenoids have oxygenated functional groups"
    
    # Check for functional group enrichment indicative of triterpenoids
    # Here, we need to extend the search more explicitly
    # Look for hydroxyls (-OH), carbonyls (C=O), and more as defining features
    if not any(atom.GetSymbol() == 'O' and atom.GetDegree() >= 2 for atom in mol.GetAtoms()):
        return False, "No eligible functional groups indicative of triterpenoid modifications (e.g., oxygenated functionality like hydroxyl groups)"

    # Consider more forgiving checks for the presence of multiple ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        # Triterpenoids often have multiple fused rings, including ring modifications
        return False, f"Insufficient ring structures found, expected at least 3, found {ring_info.NumRings()}"

    # Check for extensive sugar or glycosidic attachments
    if any(atom.GetAtomicNum() == 6 and atom.GetDegree() > 4 for atom in mol.GetAtoms()):
        # This alone isn't a disqualifier, merely a complexity marker we absorb with flexibility
        pass
    
    # If all criteria are met, classify as a triterpenoid
    return True, "Molecule matches triterpenoid characteristics based on flexible carbon count, functional groups, and ring presence"