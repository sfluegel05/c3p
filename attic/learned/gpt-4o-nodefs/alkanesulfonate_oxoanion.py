"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion contains a sulfonate group attached to an aliphatic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sulfonate group pattern
    sulfonate_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[O-]")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Missing sulfonate group S(=O)(=O)[O-]"

    # Check for aliphatic carbon connection
    aliphatic_connected = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Check if either atom is sulfur from sulfonate
        if (atom1.GetSymbol() == 'S' and atom2.GetSymbol() == 'C' and not atom2.GetIsAromatic()) or \
           (atom2.GetSymbol() == 'S' and atom1.GetSymbol() == 'C' and not atom1.GetIsAromatic()):
            # Ensure the carbon is connected to a sulfonate sulfur
            if mol.GetAtomWithIdx(atom1.GetIdx()).GetSymbol() == 'S' or mol.GetAtomWithIdx(atom2.GetIdx()).GetSymbol() == 'S':
                aliphatic_connected = True
                break
    
    if aliphatic_connected:
        return True, "Sulfonate group attached to an aliphatic structure"
    
    return False, "Sulfonate group not connected to an aliphatic structure"