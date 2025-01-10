"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer obtained by the oxidative coupling 
    of at least two units of aryl-substituted benzopyran rings or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expand pattern for flavonoid units - cover variations
    flavonoid_unit_patterns = [
        Chem.MolFromSmarts("c1cc(O)c2oc3ccccc3c(c2c1)O"), # typical flavone core
        Chem.MolFromSmarts("c1cc(O)c2c(c1)oc1ccccc1c2O"), # variational flavone structure
        Chem.MolFromSmarts("c1cc(=O)c2c(c1)c(c(c(c2O)O)O)O") # adapt for polyhydroxylated flavonoids
    ]
    match_flag = False
    flavonoid_matches_count = sum(bool(mol.HasSubstructMatch(pat)) for pat in flavonoid_unit_patterns)
    if flavonoid_matches_count < 2:
        return False, "Less than two flavonoid units found"
    
    # Improved linkage pattern to capture multiple variations
    linkage_patterns = [
        Chem.MolFromSmarts("c-C-c"),  # generic C-C linkage
        Chem.MolFromSmarts("c-c-c"),  # possible aromatic linkage
        Chem.MolFromSmarts("c-O-c"),  # oxygen bridged
        Chem.MolFromSmarts("c(=O)-c"),  # involving carbonyl linkage
    ]
    
    for linkage_pattern in linkage_patterns:
        if mol.HasSubstructMatch(linkage_pattern):
            match_flag = True
            break

    if match_flag:
        return True, "Contains flavonoid units linked by various biflavonoid patterns"
    return False, "Failed to identify expected biflavonoid linkage between flavonoid units"