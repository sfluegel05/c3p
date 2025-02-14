"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: CHEBI:51519 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify functional groups
    families = Chem.GetMolFeatureFamilies(mol)

    # Check for amino acid backbone
    has_amino_acid_backbone = False
    for family in families:
        if family.FamilyName == 'AlphaAmino':
            has_amino_acid_backbone = True
            break

    if not has_amino_acid_backbone:
        return False, "No amino acid backbone found"

    # Identify alpha carbon
    alpha_carbon = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [mol.GetAtomWithIdx(nbr.GetIdx()) for nbr in atom.GetNeighbors()]
            if any(nbr.GetSymbol() == 'N' for nbr in neighbors) and any(nbr.GetSymbol() == 'O' for nbr in neighbors):
                alpha_carbon = atom
                break

    if alpha_carbon is None:
        return False, "Unable to identify alpha carbon"

    # Check for N-hydroxy group on alpha carbon
    has_n_hydroxy = False
    for nbr in alpha_carbon.GetNeighbors():
        nbr_atom = mol.GetAtomWithIdx(nbr.GetIdx())
        if nbr_atom.GetSymbol() == 'N' and any(nnbr.GetSymbol() == 'O' for nnbr in nbr_atom.GetNeighbors()):
            has_n_hydroxy = True
            break

    if has_n_hydroxy:
        return True, "Contains amino acid backbone and N-hydroxy group on alpha carbon"
    else:
        return False, "No N-hydroxy group on alpha carbon"