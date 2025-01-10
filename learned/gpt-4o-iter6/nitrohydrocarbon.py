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

    # Verify the molecule contains a carbon backbone with hydrogens
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    hydrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1]
    if not carbon_atoms or not hydrogen_atoms:
        return False, "No sufficient carbon-hydrogen backbone found"
    
    # Check that nitro groups are substituting hydrogens on a carbon backbone
    for atom in carbon_atoms:
        if atom.GetDegree() < 4:  # A fully saturated carbon would usually have 4 bonds
            attached_nitro = any(mol.GetBondBetweenAtoms(atom.GetIdx(), n[0]) for n in nitro_matches)
            if attached_nitro:
                return True, "Molecule is a nitrohydrocarbon with nitro groups substituting hydrogens"

    return False, "Nitro groups are not part of main hydrocarbon substitution"