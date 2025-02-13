"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid contains no carbon-to-carbon multiple bonds, has a long linear carbon chain,
    and includes a terminal carboxyl group (-C(=O)O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule contains a terminal carboxyl group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Checking for a long continuous chain of carbon without any multiple bonds
    chain_length = 0
    for atom in mol.GetAtoms():
        # Count only carbon atoms in a single continuous chain
        if atom.GetAtomicNum() == 6:  # Carbon
            neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(neighbors) == 1 or len(neighbors) == 2:  # Accept single or part of a linear chain
                chain_length += 1
    
    # Ensure chain length is reasonable for fatty acids (>4 typically)
    if chain_length < 4:
        return False, f"Insufficient carbon chain length, found {chain_length} carbons"

    # Check for saturation: ensure no multiple C-C bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE,):
            begin_atom, end_atom = bond.GetBeginAtom(), bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                return False, "Contains multiple carbon-carbon bonds"

    return True, "Molecule is a saturated fatty acid with long linear carbon chain and contains a terminal carboxyl group"