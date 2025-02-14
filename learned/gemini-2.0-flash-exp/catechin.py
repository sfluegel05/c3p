"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ols, having the core structure of two benzene rings
    linked by a pyran ring and a hydroxyl group at the 3-position of the pyran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catechin, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core flavan-3-ol substructure (more flexible).
    # Allows for different bonding around the central ring, specifically,
    # the direct link to aromatic or a single bond connection with an aromatic ring.
    flavan_3ol_core = Chem.MolFromSmarts("[CH2][CH](O)[CH]([c]1[c]([c][c][c][c]1)[O]~[c]2[c][c][c][c][c]2)")
    
    if not mol.HasSubstructMatch(flavan_3ol_core):
        return False, "Core flavan-3-ol structure not found."

    # Verify the presence of two aromatic rings connected to the core structure.
    aromatic_ring = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")
    aromatic_matches = mol.GetSubstructMatches(aromatic_ring)
    if len(aromatic_matches) < 2:
        return False, "Not enough aromatic rings present for a catechin."

    # Ensure no ketone group at the 4-position (exclude flavanones)
    flavanone_core = Chem.MolFromSmarts("[CH2][C](=[O])[CH]1[c]2[c]([c][c][c][c]2)[O]1")
    if mol.HasSubstructMatch(flavanone_core):
        return False, "Molecule is a flavanone, not a catechin"
        
    
    return True, "Molecule contains the core flavan-3-ol structure and 2 aromatic rings."