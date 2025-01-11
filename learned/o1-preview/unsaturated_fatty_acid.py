"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: unsaturated fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is defined as any fatty acid containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for aliphatic chain (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a typical fatty acid"

    # Check for unsaturation: C=C or C#C bonds
    unsat_bonds = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and all(atom.GetAtomicNum() == 6 for atom in (bond.GetBeginAtom(), bond.GetEndAtom())):
            unsat_bonds = True
            break
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE and all(atom.GetAtomicNum() == 6 for atom in (bond.GetBeginAtom(), bond.GetEndAtom())):
            unsat_bonds = True
            break
    if not unsat_bonds:
        return False, "No carbon-carbon double or triple bonds found"

    # Check that the main chain is a long aliphatic chain (typical fatty acids have chain lengths >= 4)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Chain too short for fatty acid (found {num_carbons} carbons)"

    return True, "Molecule is an unsaturated fatty acid"