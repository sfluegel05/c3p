"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI: monoamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine contains at least one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for at least one aromatic ring
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
    if not aromatic_atoms:
        return False, "No aromatic ring found"

    # SMARTS pattern for aromatic -> C -> C -> N (non-aromatic nitrogen)
    # Ensures two-carbon chain between aromatic and amine group
    monoamine_pattern = Chem.MolFromSmarts("[a]-[#6]-[#6]-[N;!a]")
    matches = mol.GetSubstructMatches(monoamine_pattern)
    if not matches:
        return False, "No two-carbon chain connecting aromatic ring to amine"

    # Collect unique nitrogen atoms from matches
    amine_nitrogens = {match[-1] for match in matches}

    # Exclude nitro groups (N with adjacent O)
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if mol.HasSubstructMatch(nitro_pattern):
        return False, "Nitro group present"

    # Exclude amides (N adjacent to carbonyl)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Amide group present"

    # Check if any valid amine nitrogen remains
    valid_nitrogens = []
    for n_idx in amine_nitrogens:
        atom = mol.GetAtomWithIdx(n_idx)
        # Check if nitrogen is in a nitro/amide group (redundant but safe)
        if atom.GetAtomicNum() == 7 and not atom.IsInRing():
            valid_nitrogens.append(n_idx)

    if not valid_nitrogens:
        return False, "No valid amine groups found"

    return True, f"Found {len(valid_nitrogens)} amino group(s) connected to aromatic ring via two-carbon chain"