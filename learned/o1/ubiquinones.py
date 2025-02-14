"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define the ubiquinone core SMARTS pattern (corrected)
    # Core quinone ring with methoxy groups at positions 2 and 3, methyl at 5, and attachment at position 6
    core_smarts = "COc1cc(C)c(=O)c(OC)c1=O"
    core_mol = Chem.MolFromSmarts(core_smarts)
    if core_mol is None:
        return False, "Unable to create core pattern"

    # Find matches for the core pattern
    matches = mol.GetSubstructMatches(core_mol)
    if not matches:
        return False, "Core ubiquinone structure not found"

    # Identify the atom in the core pattern corresponding to position 6
    # In the core_smarts, the ring atoms are numbered as follows (based on SMILES parsing):
    # Atom indices in core_mol: 0 to 11 (including substituents)
    # Ring atoms indices: 2 (C1), 3, 4, 5, 6, 7
    # Position 6 corresponds to atom index 7 in core_mol

    position6_idx_in_core = 7  # Index of the ring atom at position 6

    # Check for substituent (side chain) at position 6 in the molecule
    for match in matches:
        position6_idx_in_mol = match[position6_idx_in_core]
        position6_atom = mol.GetAtomWithIdx(position6_idx_in_mol)

        # Check if the atom at position 6 has substituents outside the core structure
        core_atom_indices = set(match)
        substituents = [
            neighbor for neighbor in position6_atom.GetNeighbors()
            if neighbor.GetIdx() not in core_atom_indices
        ]

        if substituents:
            # Side chain exists; check if it is a polyprenoid chain
            side_chain_atom = substituents[0]
            side_chain_fragment = Chem.FragmentOnBonds(
                mol, [mol.GetBondBetweenAtoms(position6_idx_in_mol, side_chain_atom.GetIdx()).GetIdx()],
                addDummies=False
            )
            side_chain = Chem.GetMolFrags(side_chain_fragment, asMols=True, sanitizeFrags=True)[1]
            
            # Check if side chain is long enough (polyprenoid chains are typically long hydrocarbon chains)
            num_carbons = sum(1 for atom in side_chain.GetAtoms() if atom.GetAtomicNum() == 6)
            if num_carbons >= 5:
                return True, "Contains ubiquinone core with polyprenoid side chain at position 6"
            else:
                continue  # Side chain too short to be polyprenoid

    return False, "Core ubiquinone structure found but no suitable polyprenoid side chain at position 6"

__metadata__ = {   
    'chemical_class': {   
        'name': 'ubiquinones',
        'definition': 'Any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone; one of a group of naturally occurring homologues. The redox-active quinoid moiety usually carries a polyprenoid side chain at position 6, the number of isoprenoid units in which is species-specific. Ubiquinones are involved in the control of mitochondrial electron transport, and are also potent anti-oxidants.',
    }
}