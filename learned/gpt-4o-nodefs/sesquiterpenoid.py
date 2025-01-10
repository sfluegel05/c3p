"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is typically characterized by its 15-carbon skeleton derived from three isoprene units.
    Consideration is given to common structural motifs such as cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Approximate count of 15 carbons but allow for small variations
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12 or c_count > 18:
        return False, f"Carbon count is {c_count}, expected roughly between 12 and 18 for sesquiterpenoid"

    # Basic structural motifs can be identified using SMARTS
    # Common motifs include isoprene units, terpenic rings
    # Here we use a simplification to detect basic ring structures (e.g., cyclohexane)
    # A full set of realistic patterns needs refining.
    
    # Cyclohexane-like structure, a proxy for bicyclic sesquiterpene patterns
    cyclohexane_pattern = Chem.MolFromSmarts("C1CCCCC1")

    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No basic cyclohexane-like structure found"

    # Validate additional sesquiterpenoid-like features, like multiple rings or specific functionalities,
    # which are not straightforward with basic SMARTS
    # (Extend with more sophisticated patterns based on terpenoid biology)
    
    # This simplified model assumes basic ring presence as indicative of potential sesquiterpenoid.
    return True, "Contains structures consistent with sesquiterpenoid skeleton"