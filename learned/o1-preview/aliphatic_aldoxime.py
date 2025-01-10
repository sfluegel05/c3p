"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime derived from an aliphatic aldehyde,
    meaning it contains an aldoxime group (R-CH=N-OH) attached to an aliphatic chain,
    with no aromatic rings or cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aromatic atoms
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms or rings"
    
    # Check for any rings (acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic structures"
    
    # Define a more specific aldoxime SMARTS pattern
    # Carbon with one hydrogen, double-bonded to neutral nitrogen,
    # nitrogen single-bonded to oxygen with hydrogen
    aldoxime_pattern = Chem.MolFromSmarts("[C&H1;!R]=[N&H0;X2;!R;!+0][O&H1;X2]")
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aldoxime functional group found"
    
    # Ensure the aldoxime group is attached to an aliphatic chain (no adjacent multiple bonds)
    aliphatic_chain_pattern = Chem.MolFromSmarts("[C&H1;!R]=N[O&H1]-[C;!R;X4]")
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "Aldoxime group not attached to an aliphatic chain"
    
    return True, "Contains aldoxime group derived from an aliphatic aldehyde"