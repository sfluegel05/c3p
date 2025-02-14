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

    # Define the core aldoxime substructure using SMARTS: [C;H1,H0X3]=[N][OX2]
    # This means: A carbon atom with one explicit H or zero explicit hydrogens,
    # and 3 total connections, connected to the =N-O oxime group. Also, we explicitly say that O is attached with one bond only.
    oxime_pattern = Chem.MolFromSmarts("[C;H1,H0X3]=[N][OX2]")
    
    # Find matches for the oxime pattern
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
        return False, "Molecule does not contain the C=N-O oxime group"
    
    for match in matches:
        c_atom_idx = match[0] #Index of the carbon in the C=N
        c_atom = mol.GetAtomWithIdx(c_atom_idx)
        
        num_non_h_neighbors = 0
        num_h_neighbors = 0
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 1:
                num_non_h_neighbors += 1
            else:
                num_h_neighbors += 1
        
        #Check there is one non-hydrogen neighbor and one explicit or implicit hydrogen
        if num_non_h_neighbors == 1 and (c_atom.GetTotalNumHs() == 1):
            return True, "Molecule contains the R-CH=NOH aldoxime substructure."

    return False, "Molecule does not contain the correct R-CH=NOH aldoxime substructure."