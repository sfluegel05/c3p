"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Look for indole moiety: C1=CC2=C(NC=C2)C=C1
    indole_pattern = Chem.MolFromSmarts("c1cc2c(cc1)[nH]c2")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety found"
        
    # Check for additional nitrogen atoms (more than one in the whole molecule)
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(n_atoms) < 2:
        return False, "Not enough nitrogen atoms for alkaloid"

    # Check for complex polycyclic system
    if not Chem.GetSSSR(mol) > 5:  # Simple heuristic for complexity
        return False, "Insufficient ring complexity"

    return True, "Molecule contains an indole moiety and additional nitrogen atoms characteristic of alkaloids"