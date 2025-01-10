"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone contains an aromatic ketone with a -C=C-C(=O)- linker between two aromatic rings.

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

    # Define chalcone core pattern: aromatic ring -> C=C -> C(=O) -> aromatic ring
    chalcone_pattern = Chem.MolFromSmarts('c1ccccc1/C=C/C(=O)c2ccccc2')

    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Molecule contains the chalcone core structure"
    
    return False, "Molecule does not contain the chalcone core structure"


# Example usage:
smiles_examples = [
    "COc1cc(O)c(C(=O)\C=C\c2ccccc2)c(OC)c1",  # flavokawain B
    "CC(C)=CCc1cc(O)[C@@H](O)c(C(=O)\C=C\c2ccc(O)cc2)c1O",  # example chalcone
    "O=C(\C=C\c1ccccc1)c1ccccc1",  # trans-chalcone
    "CC(C)=CCc1c(O)ccc(C(=O)\C=C\c2ccc(O)cc2)c1O",  # isobavachalcone
]

for smiles in smiles_examples:
    is_chalcone, reason = is_chalcones(smiles)
    print(f"SMILES: {smiles} => Is Chalcone: {is_chalcone}, Reason: {reason}")