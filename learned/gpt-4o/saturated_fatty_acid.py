"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid contains no carbon-to-carbon multiple bonds and includes a carboxyl group.

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

    # Ensure the molecule contains a carboxyl group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for carbon-carbon bonds and ensure there are no double or triple bonds
    for bond in mol.GetBonds():
        begin_atom, end_atom = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE,):
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                return False, "Contains multiple carbon-carbon bonds"

    return True, "Molecule is a saturated fatty acid with no carbon-to-carbon multiple bonds and contains a carboxyl group"