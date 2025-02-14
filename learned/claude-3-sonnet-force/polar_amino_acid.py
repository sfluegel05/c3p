"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: CHEBI:33711 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# Define polar functional groups
polar_groups = ['[NH2]', '[NH3+]', '[OH]', '[OX1H]', '[SX1H]']

# Define SMARTS patterns for polar groups
polar_patterns = [Chem.MolFromSmarts(pattern) for pattern in polar_groups]

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is defined as any amino acid whose side chain is capable
    of forming one or more hydrogen bonds.

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
    aa_backbone = Chem.MolFromSmarts('N[C@H](C(=O)O)C')
    if not mol.HasSubstructMatch(aa_backbone):
        return False, "Not an amino acid"

    # Check for polar groups in side chain
    side_chain = AllChem.DeleteSubstructs(mol, aa_backbone)
    for pattern in polar_patterns:
        if side_chain.HasSubstructMatch(pattern):
            return True, "Side chain contains polar group capable of hydrogen bonding"

    return False, "Side chain does not contain polar groups"