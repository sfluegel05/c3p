"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential fatty acid 
Definition: Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
For our purposes, the molecule must be a free (non-esterified) fatty acid that is acyclic, has a single terminal
carboxylic acid group, a sufficiently long (≥10 carbon) chain, and at least 2 carbon–carbon double bonds (apart from the carbonyl).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is a (free) essential fatty acid based on its SMILES string.
    The criteria are:
      - The SMILES string must be valid.
      - The molecule must be acyclic.
      - The molecule must contain exactly one free carboxylic acid group (either protonated or ionized).
        In a free fatty acid, the carboxyl carbon must be terminal (only one carbon neighbor).
      - The molecule must contain at least 10 carbon atoms.
      - Excluding the carbonyl of the acid group, the molecule must have at least 2 carbon–carbon double bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria for an essential fatty acid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, require that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a simple unbranched fatty acid"

    # Look for a free carboxylic acid group.
    # Use two SMARTS patterns: one for the protonated acid and one for the deprotonated form.
    acid_smarts_prot = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    acid_smarts_anion = Chem.MolFromSmarts("[CX3](=O)[O-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts_prot) + mol.GetSubstructMatches(acid_smarts_anion)
    if len(acid_matches) != 1:
        return False, f"Expected one free carboxylic acid group, found {len(acid_matches)}"

    # Ensure that the acid carbon is terminal.
    # We assume that in the match the first atom is the carbonyl carbon.
    acid_carbon = mol.GetAtomWithIdx(acid_matches[0][0])
    # Count its carbon neighbors (ignoring the oxygen atoms part of the acid group).
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxylic acid group is not terminal (likely part of an ester linkage)"

    # Count the total number of carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbons)
    if num_carbons < 10:
        return False, f"Not enough carbon atoms (found {num_carbons}, need at least 10)"

    # Count the number of carbon–carbon double bonds.
    # Exclude double bonds that are part of the carboxyl group.
    cc_double_bonds = 0
    acid_carbon_idx = acid_carbon.GetIdx()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Both atoms must be carbon.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Skip the double bond if it is attached to the acid carbon (i.e. part of the COOH group)
                if a1.GetIdx() == acid_carbon_idx or a2.GetIdx() == acid_carbon_idx:
                    continue
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Not polyunsaturated enough (found {cc_double_bonds} carbon–carbon double bond(s); need at least 2)"

    return True, (f"Contains a free carboxylic acid group, {num_carbons} carbons, and {cc_double_bonds} carbon–carbon double bonds, "
                  "consistent with an essential fatty acid")