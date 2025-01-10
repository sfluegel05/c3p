"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon where one or more hydrogens have been replaced by nitro groups.
    
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
    
    # Define nitro group pattern: -[N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[NX3](=O)[O-]")

    # Check for nitrate group presence
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"

    # Verify the molecule contains carbon besides nitro group(s)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found, not a hydrocarbon base"
    
    # Count nitro groups
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    nitro_count = len(nitro_matches)
    
    # Ensure primarily a hydrocarbon with nitro substitutions
    non_carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6, 7, 8)] # N and O are part of nitro
    if len(non_carbon_atoms) > nitro_count:
        return False, "Too many non-nitro functional groups"

    return True, "Molecule is a nitrohydrocarbon with one or more nitro groups attached to a carbon backbone"