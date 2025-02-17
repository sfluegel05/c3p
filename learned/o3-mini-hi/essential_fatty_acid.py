"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential fatty acid (“any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement”)

This program checks if a given SMILES string represents a molecule that is likely
a free polyunsaturated fatty acid by using the following criteria:
  - The SMILES is valid.
  - The molecule contains exactly one carboxylic acid group ([C(=O)O]) which designates the fatty acid head.
  - The molecule has a sufficient number of carbon atoms (we require ≥10 carbons).
  - The molecule contains at least 2 carbon–carbon double bonds (i.e. excluding the carbonyl of the acid),
    consistent with being polyunsaturated.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    For our purposes, an essential fatty acid is defined as a polyunsaturated fatty acid
    displaying a free carboxylic acid group and a long, unbranched aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a polyunsaturated fatty acid satisfying our criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a carboxylic acid group.
    # The SMARTS "C(=O)[O;H,-]" should match both the neutral and ionized forms.
    acid_smarts = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    if len(acid_matches) != 1:
        return False, f"Expected one carboxylic acid group, found {len(acid_matches)}"

    # Count carbon atoms; fatty acids are long-chain molecules.
    # For our purpose, we require at least 10 carbons.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 10:
        return False, f"Not enough carbon atoms (found {len(carbon_atoms)}; need at least 10)"

    # Count the number of carbon–carbon double bonds.
    # We count bonds that are double bonds and that connect two carbon atoms.
    # The carbonyl (C=O) bond is not counted since one atom is oxygen.
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Not polyunsaturated enough (found {cc_double_bonds} carbon–carbon double bond(s); need at least 2)"

    # Optional: one might check for additional acyl chain features (e.g. unbranched, acyclic) but for simplicity
    # we consider the presence of a carboxyl group, a long carbon chain, and sufficient C=C bonds as our decision criteria.
    return True, (f"Contains a carboxylic acid group, {len(carbon_atoms)} carbons, "
                  f"and {cc_double_bonds} carbon–carbon double bonds, consistent with an essential fatty acid")