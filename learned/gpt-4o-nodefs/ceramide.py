"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import Mol

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide typically has a sphingosine backbone linked to a fatty acid via an amide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for the sphingosine backbone pattern (long chain with amine)
    sphingosine_pattern = Chem.MolFromSmarts("C[C@H](O)[C@@H](N)CO")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"

    # Search for amide linkage: [CX3](=O)[NX3]
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Verify that the amide linkage connects to a long aliphatic fatty acid chain
    connected_fatty_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    chain = rdmolops.GetLongestHomogeneousAtomPath(mol, neighbor.GetIdx(), 6)
                    if len(chain) > 10:
                        connected_fatty_acid = True
                        break
    if not connected_fatty_acid:
        return False, "Amide linkage not connected to a long fatty acid chain"

    return True, "Molecule contains a sphingosine backbone with a fatty acid linked via an amide bond"