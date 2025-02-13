"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone has a 2,3-dimethoxy-5-methylbenzoquinone core with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS for 2,3-dimethoxy-5-methylbenzoquinone core
    benzoquinone_pattern = Chem.MolFromSmarts("COc1cc(=O)c(C)c(=O)c(OC)c1")
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No 2,3-dimethoxy-5-methylbenzoquinone core found"
    
    # Find potential polyprenoid side chains with isoprene-like patterns
    long_chain_found = False
    # Check for multiple occurrences of C=C groups followed by multiple C's, characteristic of isoprene units
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = atom.GetNeighbors()
            if sum(1 for n in neighbors if n.GetSymbol() == 'C' and n.GetDegree() == 3) > 2:
                # Likely part of a polyprenoid-like structure
                long_chain_found = True
                break
    
    if not long_chain_found:
        return False, "No polyprenoid side chain found"

    return True, "Contains ubiquinone core with polyprenoid side chain"