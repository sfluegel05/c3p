"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid has no carbon-to-carbon multiple bonds, typically a long chain,
    and a terminal carboxyl group. Cyclic structures disqualify the molecule as a typical
    saturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a saturated fatty acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify terminal carboxyl group (C(=O)O)
    carboxyl_pattern = Chem.MolFromSmiles("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxyl group found"

    # Check for any carbon-carbon double or triple bonds indicating unsaturation
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                return False, "Contains carbon-carbon multiple bonds (unsaturation)"

    # Check for cyclic structures
    if Chem.GetSSSR(mol) > 0:
        return False, "Contains cyclic structure, not a typical saturated fatty acid"

    # Count carbon atoms to ensure typical fatty acid length
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count < 4:
        return False, f"Insufficient carbon chain length, found {carbon_count} carbons"

    return True, "Molecule is a saturated fatty acid with unbranched (possibly branched) carbon chain and terminal carboxyl group"