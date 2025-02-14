"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is characterized by a nitrogen atom with a +1 charge attached
    to two carbon chains, typically formed by protonation of a secondary amine (R2NH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for secondary ammonium ion
    # The nitrogen must have a +1 charge and be connected to exactly two carbon atoms, non-aromatic
    secondary_ammonium_pattern = Chem.MolFromSmarts("[NH2+;R0]([C;!R])[C;!R]")

    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        # Check each match for additional verification
        for match in mol.GetSubstructMatches(secondary_ammonium_pattern):
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            if nitrogen_atom.GetFormalCharge() == 1 and nitrogen_atom.GetDegree() == 3:
                return True, "Contains protonated secondary ammonium ion pattern"

    return False, "Does not contain a secondary ammonium ion pattern"