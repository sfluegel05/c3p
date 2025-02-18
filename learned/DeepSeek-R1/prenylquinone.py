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
    A prenylquinone has a quinone core (1,4-benzoquinone or naphthoquinone)
    substituted by at least one polyprenyl-derived side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define quinone core patterns
    benzoquinone = Chem.MolFromSmarts("[O]=C1C(=O)C=CC=C1")  # 1,4-benzoquinone
    naphthoquinone = Chem.MolFromSmarts("[O]=C1C(=O)C2=CC=CC=C2C=C1")  # 1,4-naphthoquinone

    # Check for quinone core
    has_quinone = mol.HasSubstructMatch(benzoquinone) or mol.HasSubstructMatch(naphthoquinone)
    if not has_quinone:
        return False, "No quinone core detected"

    # Define prenyl chain pattern (at least two isoprene units)
    # Matches trans-configured isoprenoid chains with methyl branches
    prenyl_pattern = Chem.MolFromSmarts("[CH2][CH](C)=[CH][CH2]")
    
    # Find all possible attachment points on quinone core
    core_matches = []
    if mol.HasSubstructMatch(benzoquinone):
        core_matches = mol.GetSubstructMatches(benzoquinone)
    else:
        core_matches = mol.GetSubstructMatches(naphthoquinone)

    # Check for prenyl substituents attached to quinone core
    core_atoms = set()
    for match in core_matches:
        core_atoms.update(match)

    # Find all prenyl chain matches
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No polyprenyl chain detected"

    # Verify at least one prenyl chain is attached to quinone core
    for pm in prenyl_matches:
        for atom_idx in pm:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if adjacent to core atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in core_atoms:
                    return True, "Quinone core with polyprenyl side chain"

    return False, "No polyprenyl chain attached to quinone core"