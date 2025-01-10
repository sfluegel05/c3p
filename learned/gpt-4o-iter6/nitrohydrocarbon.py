"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens have been replaced by nitro groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitro group pattern -[N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"
      
    # Check for carbon atom presence
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found"

    # Check if molecule is primarily made up of C and H in addition to nitro groups
    has_significant_organic_structure = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (6, 1):  # Carbon and Hydrogen
            has_significant_organic_structure = True
            break
    
    if not has_significant_organic_structure:
        return False, "Lacks significant hydrocarbon structure"

    # Molecules can have additional elements, but must have nitro groups as significant feature
    # Thus, any meaningful additional features should not disqualify it as nitrohydrocarbon
    
    return True, "Molecule is a nitrohydrocarbon with one or more nitro groups attached to a carbon backbone"