"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a monocarboxylic acid with a long aliphatic chain, containing no carbon-carbon multiple bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Check for any carbon-carbon double or triple bonds
    has_multiple_bonds = False
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                has_multiple_bonds = True
                break
    if has_multiple_bonds:
        return False, "Contains carbon-carbon multiple bonds"

    # Check if the molecule is aliphatic (no aromatic rings)
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Contains aromatic rings"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Too few carbons ({num_carbons}) for a fatty acid"

    return True, "Molecule is a saturated fatty acid"