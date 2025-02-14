"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as RCH=NOH, derived from an aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core aldoxime substructure using SMARTS: [C;H1,H0X3]=[N][O]
    # This means: A carbon atom with one explicit H or zero explicit hydrogens (implicitly one H),
    # and 3 total connections, connected to the =N-O oxime group.
    oxime_pattern = Chem.MolFromSmarts("[C;H1,H0X3]=[N][O]")
    
    # Find matches for the oxime pattern
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
      return False, "Molecule does not contain the C=N-O oxime group"

    # Check that the carbon is connected to one R group only
    for match in matches:
        c_atom_idx = match[0] #The carbon in the C=N
        c_atom = mol.GetAtomWithIdx(c_atom_idx)

        # Check the number of non-hydrogen neighbors to this carbon. If there is only one, then the pattern of RCH=NOH is met.
        num_non_h_neighbors = 0
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 1: #check if not hydrogen
                num_non_h_neighbors += 1
        if num_non_h_neighbors == 1:
            return True, "Molecule contains the R-CH=NOH aldoxime substructure."
    
    return False, "Molecule does not contain the correct R-CH=NOH aldoxime substructure."