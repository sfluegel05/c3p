"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone contains an aromatic ketone with an extended conjugated system:
    Ar-C=C-C(=O)-Ar with potential substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Chalcone pattern: aromatic ring -> C=C -> C(=O) -> aromatic ring
    # Allow for additional substitutions like OH, OCH3, etc on the rings
    chalcone_pattern_smarts = '[$([a]:[c]):!$(c1cc(C(=O)C=c2ccccc2)[cH]c1)](=C)=C(=O)[$([c]:[a]):!$(c1ccccc1C(=O)C=Cc2ccccc2)]'
    
    chalcone_pattern = Chem.MolFromSmarts(chalcone_pattern_smarts)

    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Molecule contains the chalcone core structure"
    
    return False, "Molecule does not contain the chalcone core structure"

# Example usage:
smiles_examples = [
    "COc1cc(O)c(C(=O)\C=C\c2ccccc2)c(OC)c1",  # flavokawain B
    "CC(C)=CCc1cc(O)[C@@H](O)c(C(=O)\C=C\c2ccc(O)cc2)c1O",  # example chalcone
    "O=C(\C=C\c1ccccc1)c1ccccc1",  # trans-chalcone
    "CC(C)=CCc1c(O)ccc(C(=O)\C=C\c2ccc(O)cc2)c1O",  # isobavachalcone
    "OC1=C(C(O)=C(C(=O)CCC2=CC=CC=C2)C(O)=C1)CC=C(C)C", # missed chalcone
]

for smiles in smiles_examples:
    is_chalcone, reason = is_chalcones(smiles)
    print(f"SMILES: {smiles} => Is Chalcone: {is_chalcone}, Reason: {reason}")