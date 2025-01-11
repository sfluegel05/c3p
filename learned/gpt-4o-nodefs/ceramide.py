"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem

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

    # Search for the sphingosine backbone pattern (long chain with amine and hydroxyl)
    sphingosine_pattern = Chem.MolFromSmarts("C[C@H](O)[C@@H](N)CO")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"

    # Search for amide linkage: [CX3](=O)[NX3]
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Verify that the amide linkage connects to a long aliphatic fatty acid chain
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            # Check for carbon neighbors forming a link to the amide bond
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    chain_length = 0
                    visited_atoms = set()
                    atoms_to_visit = [neighbor]
                    while atoms_to_visit and chain_length <= 20:
                        current_atom = atoms_to_visit.pop()
                        if current_atom.GetIdx() not in visited_atoms:
                            visited_atoms.add(current_atom.GetIdx())
                            if current_atom.GetSymbol() == 'C':
                                chain_length += 1
                            for n in current_atom.GetNeighbors():
                                if n.GetIdx() not in visited_atoms:
                                    atoms_to_visit.append(n)
                    if chain_length > 10:
                        return True, "Molecule contains a sphingosine backbone with a fatty acid linked via an amide bond"

    return False, "Amide linkage not connected to a long fatty acid chain"