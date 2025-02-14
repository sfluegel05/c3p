"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:33842 primary amine
A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if there is at least one nitrogen atom
    if not n_atoms:
        return False, "No nitrogen atom found"

    # Define SMARTS patterns for common functional groups
    # This is a non-exhaustive list, you can add more patterns as needed
    functional_groups = [
        Chem.MolFromSmarts("[N+]"), # Quaternary nitrogen
        Chem.MolFromSmarts("[N,O,S,P]=O"), # Nitro, sulfoxide, phosphine oxide
        Chem.MolFromSmarts("[N+](=O)[O-]"), # Nitro
        Chem.MolFromSmarts("C(=O)[O,N]"), # Carboxylic acid, amide
        Chem.MolFromSmarts("[N,O,S]C#N"), # Nitrile, isocyanate, isothiocyanate
        Chem.MolFromSmarts("C(=O)N"), # Amide
        Chem.MolFromSmarts("C(=O)O"), # Carboxylic acid, ester
        Chem.MolFromSmarts("C(=S)N"), # Thioamide
        Chem.MolFromSmarts("C(=S)S"), # Thioester
        Chem.MolFromSmarts("S(=O)(=O)"), # Sulfone
        Chem.MolFromSmarts("P(=O)"), # Phosphine oxide
        Chem.MolFromSmarts("C=N"), # Imine
    ]

    # Check for primary amine groups
    primary_amine_found = False
    for n_atom in n_atoms:
        # Check if nitrogen has at least one hydrogen and one carbon neighbor
        h_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        c_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        if h_neighbors >= 1 and c_neighbors >= 1:
            # This is a potential primary amine nitrogen
            # Check if it's part of an aromatic system
            is_aromatic = n_atom.GetIsAromatic()
            if is_aromatic:
                # Check if the aromatic nitrogen has at least one hydrogen neighbor
                aromatic_h_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1 and neighbor.GetIsAromatic())
                if aromatic_h_neighbors >= 1:
                    # This is an aromatic primary amine group
                    primary_amine_found = True
                    break
            else:
                # Check if the molecule doesn't contain any conflicting functional groups
                has_conflicting_groups = any(mol.HasSubstructMatch(group_pattern) for group_pattern in functional_groups)
                if not has_conflicting_groups:
                    primary_amine_found = True
                    break

    if primary_amine_found:
        return True, "Contains a primary amine group (-NH2) bonded to a carbon"
    else:
        return False, "No primary amine group found"