"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is a terpenoid derived from a monoterpene (C10 skeleton).
    The term includes compounds in which the C10 skeleton of the parent monoterpene
    has been rearranged or modified by the removal of one or more skeletal atoms
    (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Monoterpenes are built from two isoprene units (C5H8)
    # Attempt to find isoprene units within the molecule
    # Define isoprene unit as a 5-carbon chain with two double bonds
    isoprene_pattern = Chem.MolFromSmarts("C=C-C=C-C")
    
    # Count the number of isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)
    
    if num_isoprene_units < 2:
        # Try broader patterns to match rearranged monoterpene skeletons
        # Define monoterpene core patterns (acyclic, monocyclic, bicyclic)
        monoterpene_patterns = [
            Chem.MolFromSmarts("C1=CC=CCCC1"),        # Monocyclic
            Chem.MolFromSmarts("C1CC2CCCC1C2"),       # Bicyclic
            Chem.MolFromSmarts("C1=CC=CC=C1"),        # Benzene ring (aromatic monoterpenoids)
            Chem.MolFromSmarts("C=C(C)CCC=C(C)C")     # Acyclic monoterpene skeleton
        ]
        
        monoterpene_match = False
        for pattern in monoterpene_patterns:
            if mol.HasSubstructMatch(pattern):
                monoterpene_match = True
                break
        
        if not monoterpene_match:
            return False, "No monoterpene core structure detected"
    
    # Check for presence of oxygen atoms (typical for terpenoids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; monoterpenoids typically contain oxygen functional groups"
    
    return True, "Molecule contains monoterpene core structure with oxygen functionality"