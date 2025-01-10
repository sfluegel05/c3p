"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a benzene ring carrying exactly four chloro groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure query for benzene with four chlorines
    benzene_query = Chem.MolFromSmarts("c1ccccc1")  # Benzene ring
    chloro_query = rdqueries.AtomNumEqualsQueryAtom(17)  # Cl = 17

    # Identify benzene structures with exactly four chlorines
    for substructure in mol.GetSubstructMatches(benzene_query):
        chlorine_count = sum(1 for atom_idx in substructure if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 17)
        if chlorine_count == 4:
            return True, "Benzene ring with exactly four chloro groups found"

    return False, "No benzene ring with exactly four chloro groups found"