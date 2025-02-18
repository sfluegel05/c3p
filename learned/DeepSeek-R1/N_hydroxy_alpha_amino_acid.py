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

    # Check for carboxylic acid group (either -COOH or deprotonated form)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"

    # Verify alpha-amino acid structure (amino group on adjacent carbon)
    alpha_has_n = False
    for match in mol.GetSubstructMatches(carboxylic_acid):
        carb_carbon_idx = match[0]
        carb_carbon = mol.GetAtomWithIdx(carb_carbon_idx)
        for neighbor in carb_carbon.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                alpha_carbon = neighbor
                # Check if alpha carbon has any nitrogen neighbors
                if any(nbr.GetSymbol() == 'N' for nbr in alpha_carbon.GetNeighbors()):
                    alpha_has_n = True
                    break
            if alpha_has_n:
                break
        if alpha_has_n:
            break

    if not alpha_has_n:
        return False, "Not an alpha-amino acid"

    # Check for any N-hydroxy group (N-O single bond with H on oxygen)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                        if neighbor.GetTotalNumHs() >= 1:
                            return True, "N-hydroxy group found on an amino group"

    return False, "No N-hydroxy group found on any amino group"