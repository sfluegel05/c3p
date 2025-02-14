"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

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
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Check for any carbon-carbon double or triple bonds
    has_multiple_bonds = False
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                has_multiple_bonds = True
                break
    if has_multiple_bonds:
        return False, "Contains carbon-carbon multiple bonds"

    # Check if the molecule is aliphatic (no rings)
    ring_info = mol.GetRingInfo()
    if ring_info.IsInitialized() and ring_info.NumRings() > 0:
        return False, "Contains rings, not an aliphatic chain"

    # Count the number of carbon atoms (excluding the carboxylic acid carbon)
    carboxylic_acid_carbon_idx = carboxylic_acid_matches[0][0]
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() != carboxylic_acid_carbon_idx]
    num_carbons = len(carbon_atoms)
    if num_carbons < 3:
        return False, f"Too few carbons ({num_carbons}) for a fatty acid"

    # Check that the chain is saturated (all carbons have single bonds with carbons and hydrogens)
    for atom in carbon_atoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    return False, "Contains unsaturated bonds in aliphatic chain"

    return True, "Molecule is a saturated fatty acid"