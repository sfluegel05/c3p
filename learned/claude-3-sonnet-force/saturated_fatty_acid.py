"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:36631 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a carboxylic acid with a linear or branched carbon chain,
    containing no carbon-carbon multiple bonds in the chain.

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

    # Check for single carboxylic acid group
    carboxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and
                          sum(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()) == 1)
    if carboxyl_groups != 1:
        return False, "Must have exactly one carboxylic acid group"

    # Check for saturated carbon chain
    chain_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            neighbors = [mol.GetBondWithAtomIdx(neighbor_idx).GetBondType()
                         for neighbor_idx in atom.GetNeighbors()]
            if Chem.BondType.DOUBLE not in neighbors and Chem.BondType.TRIPLE not in neighbors:
                chain_atoms.append(atom.GetIdx())
    chain = Chem.PathToSubmol(mol, chain_atoms, atomIdxs=chain_atoms)
    if chain.GetNumAtoms() < 4:
        return False, "Carbon chain too short"

    # Check for linear or branched chain
    sssr = Chem.GetSymmSSSR(chain)
    if len(sssr) > 0:
        return False, "Cyclic structures found in carbon chain"

    return True, "Contains a linear or branched saturated carbon chain with a carboxylic acid group"