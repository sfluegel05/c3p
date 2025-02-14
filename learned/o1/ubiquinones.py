"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone is defined as a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ubiquinone core SMARTS pattern
    # Core quinone ring with methoxy groups at positions 2 and 3, methyl at 5, and side chain at 6
    core_smarts = """
    [#6]1=,:[#6]([OCH3])=,:[#6](=,:[#6]([OCH3])=,:[#6](=,:[#6]1=O)[CH3])[CH3]
    """

    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Unable to create core pattern"

    # Find matches for the core pattern
    matches = mol.GetSubstructMatches(core_pattern)
    if not matches:
        return False, "Core ubiquinone structure not found"

    # Check for polyprenoid side chain at position 6
    # Define the attachment atom (position 6) in the core pattern
    attachment_atom_idx = 5  # Index of the atom in core_pattern where side chain attaches

    for match in matches:
        # Get the atom in the molecule where the side chain should be attached
        side_chain_atom_idx = match[attachment_atom_idx]
        side_chain_atom = mol.GetAtomWithIdx(side_chain_atom_idx)

        # Get neighbors not part of the core
        core_atom_indices = set(match)
        side_chain_atoms = []

        for neighbor in side_chain_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in core_atom_indices:
                # Perform BFS to get the side chain atoms
                side_chain_frag = Chem.GetMolFrags(Chem.PathToSubmol(mol, [side_chain_atom_idx, neighbor_idx]), asMols=True)
                if side_chain_frag:
                    side_chain_mol = side_chain_frag[0]
                    # Check if the side chain is a polyprenoid
                    is_polyprenoid, reason = check_polyprenoid_side_chain(side_chain_mol)
                    if is_polyprenoid:
                        return True, "Contains ubiquinone core with polyprenoid side chain at position 6"
                    else:
                        continue

    return False, "Core ubiquinone structure found but no suitable polyprenoid side chain at position 6"

def check_polyprenoid_side_chain(side_chain_mol):
    """
    Checks if the side chain is a polyprenoid (contains isoprene units).

    Args:
        side_chain_mol (Mol): RDKit molecule of the side chain

    Returns:
        bool: True if side chain is polyprenoid, False otherwise
        str: Reason for classification
    """

    # Count the number of carbons
    num_carbons = sum(1 for atom in side_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 5:
        return False, "Side chain too short to be polyprenoid"

    # Define isoprene unit SMARTS pattern
    isoprene_smarts = "C(=C)C-C=C"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
    if isoprene_pattern is None:
        return False, "Unable to create isoprene pattern"

    # Check for isoprene units
    isoprene_matches = side_chain_mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "Side chain does not contain isoprene units"

    return True, "Side chain is polyprenoid"

__metadata__ = {   
    'chemical_class': {   
        'name': 'ubiquinones',
        'definition': 'Any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone; one of a group of naturally occurring homologues. The redox-active quinoid moiety usually carries a polyprenoid side chain at position 6, the number of isoprenoid units in which is species-specific. Ubiquinones are involved in the control of mitochondrial electron transport, and are also potent anti-oxidants.',
    }
}