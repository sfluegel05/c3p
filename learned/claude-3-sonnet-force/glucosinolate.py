"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: CHEBI:16830 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    A glucosinolate has a central C atom bonded to a glycone group via S, a sulfonated oxime group via N,
    and a side chain. The sulfate group and side chain have an anti stereochemical configuration across C=N.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycone group (e.g., glucose, rhamnose, arabinose)
    glycone_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@H]([C@@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)O")
    glycone_matches = mol.GetSubstructMatches(glycone_pattern)
    if not glycone_matches:
        return False, "No glycone group found"

    # Look for sulfate group attached to glycone via sulfur
    sulfate_pattern = Chem.MolFromSmarts("[S]C(=N[O-,S]([O-])(=O)=O)")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern, useChirality=True)
    if not sulfate_matches:
        return False, "Sulfate group not found"

    # Check if sulfate group and glycone group are connected via sulfur
    sulfate_sulfur = mol.GetAtomWithIdx(sulfate_matches[0][0])
    glycone_oxygen = mol.GetAtomWithIdx(glycone_matches[0][0])
    if not mol.GetBondBetweenAtoms(sulfate_sulfur.GetIdx(), glycone_oxygen.GetIdx()):
        return False, "Sulfate group not connected to glycone group"

    # Look for side chain attached to central carbon
    side_chain_pattern = Chem.MolFromSmarts("[C]=[N]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if not side_chain_matches:
        return False, "No side chain found"

    # Check anti stereochemistry across C=N
    central_carbon = mol.GetAtomWithIdx(side_chain_matches[0][0])
    n_atom = mol.GetAtomWithIdx(side_chain_matches[0][1])
    if central_carbon.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW or n_atom.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Incorrect stereochemistry across C=N bond"

    return True, "Contains glycone group, sulfate group, and side chain with anti stereochemistry across C=N"