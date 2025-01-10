"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate the length of the longest aliphatic carbon chain
    c_chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            dfs_stack = [(atom, 0)]
            visited = set()
            while dfs_stack:
                current_atom, length = dfs_stack.pop()
                visited.add(current_atom.GetIdx())
                length += 1
                neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and not nbr.IsInRing() and nbr.GetIdx() not in visited]
                if neighbors:
                    for nbr in neighbors:
                        dfs_stack.append((nbr, length))
                else:
                    c_chains.append(length)
    if not c_chains or max(c_chains) < 12:
        return False, f"Aliphatic chain too short ({max(c_chains) if c_chains else 0} carbons)"

    # Check for nitrogen attached to carbon (amines and amides)
    nitrogen_patterns = [
        Chem.MolFromSmarts("[CX4][NX3;$([H2]),$([H1])$([H1+]),$([H0][#6])$([H0][#1])$([H0][H0])$([H0][#1][#1])]"),  # primary amine, protonated amine
        Chem.MolFromSmarts("[CX3](=O)[NX3]")  # amide
    ]
    has_nitrogen = any(mol.HasSubstructMatch(pat) for pat in nitrogen_patterns)
    if not has_nitrogen:
        return False, "No amino or amide group attached to carbon found"

    # Check for hydroxyl groups attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups attached to carbon found"

    # Exclude peptides - check for peptide bonds (amide bonds between carbonyl and nitrogen)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) > 2:
        return False, f"Too many peptide bonds ({len(peptide_bonds)} found), possible peptide"

    # If all checks pass, classify as sphingoid
    return True, "Molecule matches sphingoid structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'sphingoid',
        'definition': 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.'
    }
}