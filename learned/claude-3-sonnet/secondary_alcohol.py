"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:35484 secondary alcohol

A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom
which has two other carbon atoms attached to it.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all alcohol groups (-OH)
    alcohol_smarts = "[OX1H]"
    alcohol_groups = mol.GetSubstructMatches(Chem.MolFromSmarts(alcohol_smarts))

    for alcohol_idx in alcohol_groups:
        alcohol_atom = mol.GetAtomWithIdx(alcohol_idx)
        carbon_neighbors = [nbr.GetIdx() for nbr in alcohol_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]

        # Check if carbon has exactly 2 other carbon neighbors (saturated secondary carbon)
        if len(carbon_neighbors) == 2:
            carbon1 = mol.GetAtomWithIdx(carbon_neighbors[0])
            carbon2 = mol.GetAtomWithIdx(carbon_neighbors[1])

            # Check if both carbons are saturated
            if sum(bond.GetBondType() == Chem.BondType.SINGLE for bond in carbon1.GetBonds()) == carbon1.GetDegree() and \
               sum(bond.GetBondType() == Chem.BondType.SINGLE for bond in carbon2.GetBonds()) == carbon2.GetDegree():
                return True, "Contains a secondary alcohol group (-OH attached to saturated carbon with 2 other carbon neighbors)"

    return False, "No secondary alcohol group found"