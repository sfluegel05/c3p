"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is defined by having a steroid backbone with a beta-configured hydroxy group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General steroid backbone pattern - four rings A/B/C/D with flexibility in stereo and connectivity
    steroid_patterns = [
        Chem.MolFromSmarts("[#6]12CC[C@H]3[C@H](CCC4=C(C=CC=C4)C3)CCC1C2"),  # typical steroid pattern
        Chem.MolFromSmarts("[#6]12CC[C@@H]3[C@H](CCC4=CCCCC34)CCC1C2")  # variant with slight deviations
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No suitable steroid backbone found"

    # Look for the 17beta-hydroxy group configuration specifically
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4:  # Look for quaternary carbon
            neighbors = atom.GetNeighbors()
            for nei in neighbors:
                if nei.GetAtomicNum() == 8 and nei.GetChiralTag() in [
                    Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW
                ]:  # Check for hydroxy with stereochemistry
                    # This orientation check assumes C17 has a distinct stereochemical marker
                    if "beta" in nei.GetProp('_CIPCode', '').lower():
                        return True, "17beta hydroxy group confirmed with stereochemistry"

    return False, "No 17beta-hydroxy steroid configuration detected"