"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:50747 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phenol pattern (benzene ring with exactly one hydroxyl group)
    phenol_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1[OH]")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    
    # Ensure there is exactly one hydroxyl group on the benzene ring
    if len(phenol_matches) != 1:
        return False, f"Found {len(phenol_matches)} hydroxyl groups on benzene ring, need exactly 1"

    # Look for nitro groups (-NO2)
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"

    # Ensure the nitro group is attached to the benzene ring
    nitro_atom = nitro_matches[0][0]  # Nitrogen atom of the nitro group
    benzene_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRing() and atom.GetAtomicNum() == 6]
    
    # Check if the nitro group is bonded to a carbon in the benzene ring
    nitro_bonded_to_benzene = False
    for neighbor in mol.GetAtomWithIdx(nitro_atom).GetNeighbors():
        if neighbor.GetIdx() in benzene_atoms:
            nitro_bonded_to_benzene = True
            break

    if not nitro_bonded_to_benzene:
        return False, "Nitro group is not attached to the benzene ring"

    return True, "Contains a phenol group with exactly one nitro substituent on the benzene ring"