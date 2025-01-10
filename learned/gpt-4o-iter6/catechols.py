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

    # Define a more flexible catechol SMARTS pattern
    catechol_pattern = Chem.MolFromSmarts("c1(O)cccc(O)c1")  # Simpler representation focusing on ortho substitution
    
    # Check for catechol substructure
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains a catechol moiety (o-diphenol component)"
    else:
        # Try a more generalized path for other aromatic layers if flexibility fails
        aromatic_ring = Chem.MolFromSmarts("c1ccccc1")  # General aromatic ring
        oxy_pattern = Chem.MolFromSmarts("[OH2]")  # More flexible hydroxyl to match oxygens more easily
        n_aromatic_rings = len(mol.GetSubstructMatches(aromatic_ring))
        
        # Criterion: At least one aromatic ring and ortho adjacent hydroxyls
        hydroxyl_matches = mol.GetSubstructMatches(oxy_pattern)
        for idx, hydroxyl in enumerate(hydroxyl_matches[:-1]):
            if mol.GetBondBetweenAtoms(hydroxyl[0], hydroxyl_matches[idx + 1][0]) is not None:
                if n_aromatic_rings > 0:
                    return True, "Contains overlapping hydroxyls on an aromatic ring (catechol moiety)"
                
        return False, "No catechol moiety found"

# Examples and usage
examples = [
    "O[C@H]([C@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C(O)=O)C(O)=O",
    "Oc1cc(O)cc(O)c1"
]
for example in examples:
    result, reason = is_catechols(example)
    print(f"SMILES: {example} -> Is Catechol: {result}, Reason: {reason}")