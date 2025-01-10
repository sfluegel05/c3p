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

    # Check for the presence of nitro groups
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"

    # Verify the molecule contains carbon besides nitro group(s)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found, not a hydrocarbon base"
    
    # Count the nitro groups and ensure they are part of hydrocarbon substitution
    nitro_count = len(nitro_matches)
    
    # Check for non-hydrocarbon atoms excluding the nitro groups (exclude C, H, N, O)
    non_hydrocarbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6, 7, 8)]
    if len(non_hydrocarbon_atoms) > 0:
        return False, "Other functional groups overshadow hydrocarbon nature"

    # Ensure a significant portion of the molecule is hydrocarbon with nitro groups attached
    if nitro_count > 0 and len(carbon_atoms) >= nitro_count:
        return True, "Molecule is a nitrohydrocarbon with nitro groups substituting one or more hydrogens on a carbon skeleton"

    return False, "Insufficient hydrocarbon structure or imbalance with nitro groups"