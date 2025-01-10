"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid typically contains a steroid backbone, a ketone group on the third carbon, 
    and a double bond between carbons 4 and 5.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broadened steroid backbone pattern to include flexible stereochemistry
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C=O)CCCC3C4CCC(=O)C=C4C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # 3-oxo group typically on the A ring (strict yet accounts for substitution locations)
    oxo_pattern = Chem.MolFromSmarts("[CR1]=O")
    oxo_matches = [match for match in mol.GetSubstructMatches(oxo_pattern) if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match)]
    if not oxo_matches:
        return False, "Missing 3-oxo group"

    # The Delta(4) double bond pattern may be more specifically tuned to general knowledge
    delta4_pattern = Chem.MolFromSmarts("C=C[C@H]([C@H]1CCC2=CC(=O)CC[C@]2(C)C1)[C@H]1CCC[C@]1(O)C")
    if not mol.HasSubstructMatch(delta4_pattern):
        return False, "Missing Delta(4) double bond"

    return True, "Contains the structural features of a 3-oxo-Delta(4) steroid"