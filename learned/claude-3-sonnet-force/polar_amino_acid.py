"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: CHEBI:26976 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# Define SMARTS patterns for amino acid backbone and polar groups
AA_BACKBONE = Chem.MolFromSmarts("N[C@@H](C(=O)O)C")  # L-amino acid backbone
POLAR_GROUPS = [
    Chem.MolFromSmarts("[OH1]"),  # Hydroxyl group
    Chem.MolFromSmarts("[NH2]"),  # Primary amine
    Chem.MolFromSmarts("[NH1]"),  # Secondary amine
    Chem.MolFromSmarts("[SH1]"),  # Thiol group
    Chem.MolFromSmarts("[NH]=C([NH2])[NH2+]"),  # Guanidino group
    Chem.MolFromSmarts("c1cnc[nH]1"),  # Imidazole ring
]

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is an amino acid whose side chain is capable of forming
    one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid backbone
    if not mol.HasSubstructMatch(AA_BACKBONE):
        return False, "No amino acid backbone found"

    # Check for polar groups in the side chain
    side_chain = Chem.DeleteSubstructs(mol, AA_BACKBONE)
    has_polar_group = any(side_chain.HasSubstructMatch(pattern) for pattern in POLAR_GROUPS)

    if has_polar_group:
        return True, "Contains amino acid backbone and polar side chain"
    else:
        return False, "Non-polar side chain"

    # Handling zwitterionic forms
    # ...

    # Additional checks or filtering
    # ...