"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol where the hydroxyl group is attached to a non-aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Check that the alcohol is attached to an aliphatic carbon
    for match in mol.GetSubstructMatches(alcohol_pattern):
        for atom_idx in match:
           
            # Get the atom connected to the oxygen
            o_atom = mol.GetAtomWithIdx(atom_idx)
            neighbor_atoms = o_atom.GetNeighbors()
            if len(neighbor_atoms) != 1:
                return False, "Oxygen not bound to a single carbon"

            carbon_atom = neighbor_atoms[0]
            if carbon_atom.GetAtomicNum() != 6:
                return False, "Alcohol not bound to a carbon"

            if carbon_atom.GetIsAromatic():
                 return False, "Alcohol group is bound to an aromatic carbon"
            
    
    return True, "Molecule contains at least one alcohol group attached to a non-aromatic carbon."