"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: CHEBI:38154 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is defined as any amino acid whose side chain is capable of forming one or more hydrogen bonds.

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

    # Check if the molecule contains an amino group (-NH2) and a carboxyl group (-COOH)
    amino_pattern = Chem.MolFromSmarts("N")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(amino_pattern) or not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Not an amino acid (missing amino or carboxyl group)"

    # Check for polar side chains
    polar_patterns = [
        Chem.MolFromSmarts("[OH1]"),  # Hydroxyl group
        Chem.MolFromSmarts("[SH1]"),  # Thiol group
        Chem.MolFromSmarts("[NH2]"),  # Primary amine group
        Chem.MolFromSmarts("[NH1]"),  # Secondary amine group
        Chem.MolFromSmarts("[NH0+]"),  # Quaternary amine group
        Chem.MolFromSmarts("[OX2H1]"),  # Carboxyl group
        Chem.MolFromSmarts("[cX3]1[cX3H][nX3H][cX3H][cX3H][cX3H]1"),  # Imidazole ring (histidine)
    ]
    is_polar = any(mol.HasSubstructMatch(pattern) for pattern in polar_patterns)

    if is_polar:
        return True, "Amino acid with a polar side chain capable of forming hydrogen bonds"
    else:
        return False, "Non-polar amino acid"