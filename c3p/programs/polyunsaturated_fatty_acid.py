"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated Fatty Acid
Definition: Any fatty acid containing more than one double bond.
These molecules contain a carboxylic acid group and are reported to have cardioprotective effects.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is defined here as a fatty acid (must have a carboxylic acid group)
    that contains more than one non-aromatic carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule qualifies as a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a carboxylic acid group. We look for a pattern: C(=O)[O;H]
    carboxylic_acid_smarts = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(carboxylic_acid_smarts):
        return False, "No carboxylic acid group found; not a fatty acid"

    # Count the number of non-aromatic carbon-carbon double bonds.
    double_bond_count = 0
    for bond in mol.GetBonds():
        # Check if the bond is a double bond.
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Exclude aromatic double bonds and confirm both atoms are carbons.
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if not bond.GetIsAromatic() and atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                double_bond_count += 1

    if double_bond_count <= 1:
        return False, (f"Found {double_bond_count} carbon-carbon double bond(s); "
                       "need more than one to qualify as polyunsaturated")

    # Optional: check if the molecule has a sufficient number of carbons to be considered a long-chain fatty acid.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 8:
        return False, "The carbon chain is too short to be a typical fatty acid"

    return True, (f"Contains a carboxylic acid group and {double_bond_count} "
                  "non-aromatic carbon-carbon double bonds; qualifies as a polyunsaturated fatty acid")