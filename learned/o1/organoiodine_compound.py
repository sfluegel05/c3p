"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: organoiodine compound
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine (C-I) bond where iodine
    is directly bonded to carbon and not bonded to any heteroatoms other than carbon or hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for carbon directly bonded to iodine
    ci_pattern = Chem.MolFromSmarts("[#6]-[#53]")
    matches = mol.GetSubstructMatches(ci_pattern)

    # Iterate over all matches of the pattern
    for match in matches:
        c_idx, i_idx = match
        iodine_atom = mol.GetAtomWithIdx(i_idx)
        carbon_atom = mol.GetAtomWithIdx(c_idx)

        # Check if iodine is only bonded to carbon and hydrogen
        iodine_neighbors = [nbr.GetAtomicNum() for nbr in iodine_atom.GetNeighbors()]
        # Remove duplicates and exclude the carbon it is bonded to
        iodine_neighbors = [num for num in iodine_neighbors if num != 6]

        # If iodine is only bonded to carbon and hydrogen (atomic number 1)
        if all(num in [1] for num in iodine_neighbors):
            return True, "Contains carbon-iodine bond with iodine in standard oxidation state (-1)"

    return False, "No suitable carbon-iodine bonds found"