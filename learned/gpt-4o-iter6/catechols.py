"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol contains an o-diphenol component, which is a benzene ring with two adjacent hydroxyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol SMARTS pattern
    # Adjust pattern to better capture variants and adjacent aromatic hydroxyls
    catechol_patterns = [
        Chem.MolFromSmarts("c1c(O)c(O)ccc1"), # Basic catechol pattern
        Chem.MolFromSmarts("c1cc(O)c(O)cc1"), # Alternate flexible pattern
        Chem.MolFromSmarts("c1(O)c(O)cccc1")  # Consider minor ring position variations
    ]

    for pattern in catechol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a catechol moiety (o-diphenol component)"
    
    return False, "No catechol moiety found"

# Test with various SMILES strings
test_smiles = {
    "O[C@H]([C@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C(O)=O)C(O)=O": "(2S,3R)-trans-caftaric acid",
    "Oc1cc(O)cc(O)c1": "Simple catechol",
    "C=1(C=CC(=C(C1)O)O)/C=C/C(OCC)=O": "Ethyl trans-caffeate"
}

for smiles, name in test_smiles.items():
    result, reason = is_catechols(smiles)
    print(f"SMILES: {smiles} ({name}) -> Is Catechol: {result}, Reason: {reason}")