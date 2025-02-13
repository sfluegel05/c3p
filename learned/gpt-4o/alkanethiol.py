"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as having a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for thiol group attached specifically to an alkyl chain
    # Sulfur atom is connected to hydrogen and an sp3 carbon (carbon without multiple bonds)
    thiol_pattern = Chem.MolFromSmarts("[SX2H1][CH2]")  # SH must be attached to sp3 carbon

    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol group (-SH) attached to an alkyl group found"

    # Extra exclusion for complex cases or attached functionalities that are not aliphatic
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            connected_atoms = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            # Check if sulfur only bonds with allowable elements for typical alkanethiols
            # Avoid connections with typical peptide/amide forming elements like nitrogen (7)
            if 7 in connected_atoms:
                return False, "Sulfur atom bonded with non-alkyl group components, possibly involving peptide or amide bonds"

    return True, "Contains a sulfanyl group, -SH, attached to an alkyl chain"