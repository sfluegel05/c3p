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
    
    # Define catechol SMARTS pattern to detect o-diphenol groups on aromatic rings
    catechol_pattern = Chem.MolFromSmarts("c1c(O)ccc(O)c1")  # Captures ortho-dihydroxy on aromatic
    
    # Check for catechol substructure
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains a catechol moiety (o-diphenol component)"
    
    # If the specific pattern fails, try a flexible approach on adjacency
    aromatic_atom_indices = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    hydroxyl_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and mol.GetAtomWithIdx(atom.GetIdx()).GetNeighbors()[0].GetIsAromatic()]

    # Check adjacency of hydroxyl groups on aromatic systems
    for i, oxygen1 in enumerate(hydroxyl_indices):
        for oxygen2 in hydroxyl_indices[i+1:]:
            if mol.GetBondBetweenAtoms(oxygen1, oxygen2):
                return True, "Contains ortho-hydroxyl groups on an aromatic ring"
    
    return False, "No catechol moiety found"

# Test with different SMILES strings
test_smiles = [
    "O[C@H]([C@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C(O)=O)C(O)=O",
    "Oc1cc(O)cc(O)c1",
    "C=1(C=CC(=C(C1)O)O)/C=C/C(OCC)=O"
]

for smiles in test_smiles:
    result, reason = is_catechols(smiles)
    print(f"SMILES: {smiles} -> Is Catechol: {result}, Reason: {reason}")