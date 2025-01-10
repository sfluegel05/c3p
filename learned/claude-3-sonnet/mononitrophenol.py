"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:51088 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at unspecified position.

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

    # More general phenol pattern that allows substitutions
    phenol_pattern = Chem.MolFromSmarts("c1([OH])c([*,H])c([*,H])c([*,H])c([*,H])c1([*,H])")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol group found"

    # Check for nitro groups - both resonance forms
    nitro_patterns = [
        Chem.MolFromSmarts("[N+](=[O-])=O"),  # Charged form
        Chem.MolFromSmarts("N(=O)=O"),         # Uncharged form
        Chem.MolFromSmarts("[N+]([O-])=O")     # Alternative charged form
    ]
    
    total_nitro = 0
    for pattern in nitro_patterns:
        if pattern is not None:  # Ensure pattern is valid
            matches = len(mol.GetSubstructMatches(pattern))
            total_nitro += matches

    if total_nitro == 0:
        return False, "No nitro group found"
    elif total_nitro > 1:
        return False, f"Found {total_nitro} nitro groups, must have exactly one"

    # Check if nitro group is attached to the phenol ring
    # More flexible pattern that allows for different positions
    nitrophenol_patterns = [
        Chem.MolFromSmarts("c1([OH])c([*,H])c([*,H])c([*,H])c([*,H])c1([N+]([O-])=O)"),
        Chem.MolFromSmarts("c1([OH])c([N+]([O-])=O)c([*,H])c([*,H])c([*,H])c1([*,H])"),
        Chem.MolFromSmarts("c1([OH])c([*,H])c([N+]([O-])=O)c([*,H])c([*,H])c1([*,H])"),
        Chem.MolFromSmarts("c1([OH])c([*,H])c([*,H])c([N+]([O-])=O)c([*,H])c1([*,H])"),
        # Add uncharged form patterns
        Chem.MolFromSmarts("c1([OH])c([*,H])c([*,H])c([*,H])c([*,H])c1(N(=O)=O)"),
        Chem.MolFromSmarts("c1([OH])c(N(=O)=O)c([*,H])c([*,H])c([*,H])c1([*,H])"),
        Chem.MolFromSmarts("c1([OH])c([*,H])c(N(=O)=O)c([*,H])c([*,H])c1([*,H])"),
        Chem.MolFromSmarts("c1([OH])c([*,H])c([*,H])c(N(=O)=O)c([*,H])c1([*,H])")
    ]

    for pattern in nitrophenol_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, "Contains phenol ring with single nitro substituent"

    return False, "Nitro group not properly attached to phenol ring"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51088',
        'name': 'mononitrophenol',
        'definition': 'A nitrophenol that is phenol carrying a single nitro substituent at unspecified position.',
        'parents': ['CHEBI:33622']
    }
}