"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid typically contains an indole moiety and multiple nitrogen atoms in a complex ring system.

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
    
    # Look for indole moiety: Structure variability included (atoms may vary)
    indole_pattern = Chem.MolFromSmarts('c1ccc2[nH]c3c(cccc3)c2c1')
    alternate_indole_pattern = Chem.MolFromSmarts('c1nccc2ccccc12')  # Alternative indole representation
    if not mol.HasSubstructMatch(indole_pattern) and not mol.HasSubstructMatch(alternate_indole_pattern):
        return False, "No indole moiety found"
        
    # Check for additional nitrogen atoms (typically more than one in the whole molecule)
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(n_atoms) < 2:
        return False, "Not enough nitrogen atoms for alkaloid"

    # Check for complex polycyclic system, using RDKit complexity descriptor
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:  # Assume indole alkaloids to have at least 3 rings
        return False, "Insufficient ring complexity"

    # Count stereocenters as another heuristic for complexity (indole alkaloids often have them)
    n_stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    if len(n_stereocenters) < 1:  # Often these have at least one stereocenter
        return False, "Insufficient stereochemistry complexity"

    return True, "Molecule contains an indole moiety, additional nitrogen, and complex ring structure"