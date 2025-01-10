"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid typically contains an indole moiety along with nitrogen atoms
    integrated into complex polycyclic structures with stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended Indole detection with varied SMARTS patterns
    indole_patterns = [
        Chem.MolFromSmarts('c1ccc2[nH]c3c(cccc3)c2c1'),  # Standard indole
        Chem.MolFromSmarts('c1nccc2ccccc12'),           # Alternate indole representation
        Chem.MolFromSmarts('c1cnc2ccccc12'),            # Nitrogen in another position
        Chem.MolFromSmarts('n1c2ccccc2cc1')             # More variations
    ]
    
    # Check if any indole pattern matches
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No indole moiety found"
        
    # Check for additional nitrogen atoms (typically more than one in a polycyclic system)
    n_atoms_ring = [atom for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() == 7 and atom.IsInRing()]
    if len(n_atoms_ring) < 2:
        return False, "Not enough nitrogen atoms in ring system"

    # Check for complex polycyclic system
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, "Insufficient ring complexity"

    # Check for stereochemistry to confirm complexity
    n_stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    if len(n_stereocenters) < 1:
        return False, "Insufficient stereochemistry complexity"

    return True, "Molecule contains an indole moiety, additional nitrogen in rings, and complex ring structure"