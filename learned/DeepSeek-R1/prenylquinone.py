"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI: prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    Requires a 1,4-benzoquinone or 1,4-naphthoquinone core with at least one
    polyprenyl-derived side chain containing multiple isoprene units.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Improved quinone core patterns with tolerance for substitutions
    # 1,4-benzoquinone: O=C1C(=O)C=CC=C1 with possible substitutions
    benzoquinone = Chem.MolFromSmarts("[O]=C1C(=O)C=CC=C1")
    # 1,4-naphthoquinone: O=C1C(=O)C2=CC=CC=C2C=C1 (naphthalene-1,4-dione)
    naphthoquinone = Chem.MolFromSmarts("[O]=C1C(=O)C2=CC=CC=C2C=C1")

    # Check for quinone core existence
    has_quinone = mol.HasSubstructMatch(benzoquinone) or mol.HasSubstructMatch(naphthoquinone)
    if not has_quinone:
        return False, "No quinone core detected"

    # Enhanced prenyl chain pattern: look for characteristic branching pattern
    # Matches at least 3 consecutive isoprene units (C-C(=C)-C-C(=C)-C-C(=C))
    prenyl_pattern = Chem.MolFromSmarts(
        "[CH3]-[CH2]-[CH]([CH3])-"
        "[CH2]-[CH]([CH3])-"
        "[CH2]-[CH]([CH3])"
    )

    # Find all possible attachment points on quinone core
    core_matches = mol.GetSubstructMatches(benzoquinone) or mol.GetSubstructMatches(naphthoquinone)

    # Check all possible prenyl chain attachments (direct or via oxygen)
    for core_match in core_matches:
        core_atoms = set(core_match)
        # Look for chains connected to core via single bond (including through oxygen)
        for atom_idx in core_atoms:
            core_atom = mol.GetAtomWithIdx(atom_idx)
            # Check neighbors for potential chain attachments
            for neighbor in core_atom.GetNeighbors():
                # Allow connection through oxygen (e.g., methoxy-linked chains)
                if neighbor.GetAtomicNum() in [6,8]:
                    # Traverse the chain looking for prenyl pattern
                    chain_root = neighbor if neighbor.GetAtomicNum() == 6 else None
                    if chain_root:
                        # Check if any part of this branch matches prenyl pattern
                        if mol.HasSubstructMatch(prenyl_pattern):
                            # Verify chain length (at least 15 carbons for 3 isoprene units)
                            chain_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
                            if chain_carbons >= 15:
                                return True, "Quinone core with polyprenyl side chain"

    # Alternative check using substructure matching with flexible attachment
    if mol.HasSubstructMatch(prenyl_pattern):
        # Verify chain is actually attached to quinone core
        combined_pattern = Chem.MolFromSmarts(
            "([O]=C1C(=O)C=CC=C1.[CH3]-[CH2]-[CH]([CH3])-[CH2]-[CH]([CH3])-[CH2]-[CH]([CH3])]) |" +
            "([O]=C1C(=O)C2=CC=CC=C2C=C1.[CH3]-[CH2]-[CH]([CH3])-[CH2]-[CH]([CH3])-[CH2]-[CH]([CH3])])"
        )
        if mol.HasSubstructMatch(combined_pattern):
            return True, "Quinone core with polyprenyl side chain"

    return False, "No polyprenyl chain attached to quinone core"