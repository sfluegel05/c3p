"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    An N-hydroxy-alpha-amino-acid is an alpha-amino acid with at least one N-hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for carboxylic acid group (either -COOH or deprotonated)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"

    # Find alpha carbon (adjacent to carboxylic acid's carbonyl carbon)
    alpha_carbon = None
    for match in mol.GetSubstructMatches(carboxylic_acid):
        carbonyl_carbon = match[0]
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                alpha_carbon = neighbor
                break
        if alpha_carbon:
            break
    if not alpha_carbon:
        return False, "No alpha carbon found"

    # Check if alpha carbon has a nitrogen attached
    nitrogen = None
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetSymbol() == 'N':
            nitrogen = nbr
            break
    if not nitrogen:
        return False, "No nitrogen attached to alpha carbon"

    # Check for N-hydroxy group (N-O single bond)
    has_n_oxide = False
    for nbr in nitrogen.GetNeighbors():
        if nbr.GetSymbol() == 'O':
            bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                has_n_oxide = True
                break
    if not has_n_oxide:
        return False, "No N-hydroxy group found"

    # Exclude amides/hydroxamates (nitrogen connected to carbonyl carbon)
    for nbr in nitrogen.GetNeighbors():
        if nbr.GetSymbol() == 'C':
            for bond in nbr.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
                    if any(atom.GetSymbol() == 'O' for atom in atoms):
                        return False, "Nitrogen is part of amide/hydroxamate group"

    return True, "N-hydroxy group found on amino group attached to alpha carbon"