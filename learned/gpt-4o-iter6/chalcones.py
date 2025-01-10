"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone contains an aromatic ketone with an extended conjugated system:
    Ar-CH=CH-C(=O)-Ar with potential substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for chalcone core structure with flexibility
    # Pattern: aromatic ring -> single bond or double bond -> C(=O) -> aromatic ring
    chalcone_pattern_smarts = 'c1ccccc1[C]=[C]-[C](=O)-c2ccccc2'
    
    # Create a molecule from the SMARTS pattern
    chalcone_pattern = Chem.MolFromSmarts(chalcone_pattern_smarts)
    if chalcone_pattern is None:
        return (None, "Unable to construct chalcone pattern.")

    # Check if the molecule matches the chalcone SMARTS pattern
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Molecule contains the chalcone core structure"
    
    return False, "Molecule does not contain the chalcone core structure"

# Example usage
smiles_examples = [
    "COc1cc(O)c(C(=O)\\C=C\\c2ccccc2)c(OC)c1",  # flavokawain B
    "OC=1C(C(C)(C)C=C)=C(O)C=CC1C(=O)/C=C/C2=CC=CC=C2",  # 2',4'-Dihydroxy-3'-(1,1-dimethyl-2-propenyl)chalcone
    "O=C(\\C=C\\c1ccccc1)c1ccccc1",  # trans-chalcone
    "OC(CC1=CC(CC=C(C)C)=C(O)C=C1)C(=O)C=2C=C(CC=C(C)C)C(O)=CC2O",  # (R)-Kanzonol Y
]

for smiles in smiles_examples:
    is_chalcone, reason = is_chalcones(smiles)
    print(f"SMILES: {smiles} => Is Chalcone: {is_chalcone}, Reason: {reason}")