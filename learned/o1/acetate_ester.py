"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdchem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for ester groups: carbonyl carbon connected to carbonyl oxygen (=O) and ester oxygen (-O)
    ester_pattern = Chem.MolFromSmarts('C(=O)O[!$([O-])]')
    if ester_pattern is None:
        return False, "Invalid SMARTS pattern for ester group"

    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    acetate_ester_count = 0

    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        carbonyl_o_idx = match[1]  # Carbonyl oxygen index
        ester_o_idx = match[2]     # Ester oxygen index

        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Get neighbors of carbonyl carbon excluding the oxygens
        neighbor_atoms = [
            nbr for nbr in carbonyl_c.GetNeighbors()
            if nbr.GetIdx() not in (carbonyl_o_idx, ester_o_idx)
        ]

        # Check if the carbonyl carbon is connected to a methyl group
        is_acetate = False
        for neighbor_atom in neighbor_atoms:
            # The neighbor should be a carbon with exactly three hydrogens (methyl group)
            if (neighbor_atom.GetAtomicNum() == 6 and
                neighbor_atom.GetTotalDegree() == 1 and
                neighbor_atom.GetTotalNumHs(includeNeighbors=True) == 3):
                is_acetate = True
                break

        if is_acetate:
            acetate_ester_count += 1

    if acetate_ester_count > 0:
        return True, f"Contains {acetate_ester_count} acetate ester group(s)"
    else:
        return False, "No acetate ester groups found"