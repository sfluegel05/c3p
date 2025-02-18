"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:46245 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone
    with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure pattern: 2,3-dimethoxy-5-methylbenzoquinone
    # with a substituent (side chain) at position 6
    core_pattern = Chem.MolFromSmarts("[#6]1(=O)[#6](-O-[#6])(-[#8]-[#6])[#6](=O)[#6](-[#6])[#6](-[!#1])[#6]1")
    if not core_pattern:
        return False, "Failed to parse core SMARTS"

    # Check for core structure match
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Core benzoquinone structure not found"

    # Verify methoxy groups (2 and 3), methyl (5), and side chain (6)
    # Additional checks to confirm substituents
    methoxy_count = 0
    methyl_found = False
    side_chain_found = False

    # Iterate over atoms in the core match
    for match in core_matches:
        ring_atoms = match
        # Assuming the core pattern's atoms are ordered as per the SMARTS
        # Adjust indices based on SMARTS structure
        # Positions 2 and 3: methoxy groups (atoms 1 and 2 in the SMARTS)
        # Position 5: methyl (atom 4 in the SMARTS)
        # Position 6: side chain (atom 5 in the SMARTS)
        # Note: SMARTS atom indices may vary; this requires testing
        # This is a simplified check and may need adjustment
        methoxy_atom1 = mol.GetAtomWithIdx(ring_atoms[1])
        methoxy_atom2 = mol.GetAtomWithIdx(ring_atoms[2])
        methyl_atom = mol.GetAtomWithIdx(ring_atoms[4])
        side_chain_atom = mol.GetAtomWithIdx(ring_atoms[5])

        # Check methoxy groups
        for atom in [methoxy_atom1, methoxy_atom2]:
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetSymbol() == 'O':
                        o_neighbors = [n.GetSymbol() for n in neighbor.GetNeighbors()]
                        if 'C' in o_neighbors and len(o_neighbors) == 1:
                            methoxy_count += 1

        # Check methyl group
        methyl_neighbors = [n.GetSymbol() for n in methyl_atom.GetNeighbors() if n.GetSymbol() == 'C']
        if len(methyl_neighbors) >= 1:
            methyl_found = True

        # Check side chain (any non-hydrogen substituent)
        if side_chain_atom.GetDegree() > 1:
            side_chain_found = True

    if methoxy_count >= 2 and methyl_found and side_chain_found:
        return True, "Matches ubiquinone core structure with methoxy, methyl, and side chain"
    else:
        return False, "Missing required substituents on core structure"