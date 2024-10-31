from rdkit import Chem
from rdkit.Chem import AllChem

def is_methylindole(smiles: str):
    """
    Determines if a molecule is a methylindole (indole with one or more methyl substituents).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a methylindole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS patterns for indole core and common variants
    indole_patterns = [
        # Basic indole cores
        Chem.MolFromSmarts("c12[nH]ccc1cccc2"),  # Unsubstituted NH indole
        Chem.MolFromSmarts("c12n(*)ccc1cccc2"),   # N-substituted indole
        Chem.MolFromSmarts("c12[nH]c(*)c1cccc2"), # 2-substituted indole
        Chem.MolFromSmarts("c12[nH]cc(*)c1cccc2"), # 3-substituted indole
        
        # Oxindole variants
        Chem.MolFromSmarts("O=C1Nc2ccccc2C1"), # Basic oxindole
        Chem.MolFromSmarts("O=C1Nc2ccccc2C1(*)"), # 3-substituted oxindole
    ]
    
    # Check for indole/oxindole core
    has_core = False
    core_atoms = set()
    for pattern in indole_patterns:
        if mol.HasSubstructMatch(pattern):
            has_core = True
            core_atoms = set(mol.GetSubstructMatch(pattern))
            break
            
    if not has_core:
        return False, "No indole core found"

    # SMARTS pattern for methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    
    # Find all methyl groups
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if not methyl_matches:
        return False, "No methyl substituents found"
        
    # Count methyl groups attached to indole/oxindole core
    methyl_count = 0
    for methyl in methyl_matches:
        methyl_carbon = mol.GetAtomWithIdx(methyl[0])
        for neighbor in methyl_carbon.GetNeighbors():
            if neighbor.GetIdx() in core_atoms or any(n.GetIdx() in core_atoms for n in neighbor.GetNeighbors()):
                methyl_count += 1
                break
                
    if methyl_count == 0:
        return False, "No methyl groups connected to indole core"
        
    return True, f"Methylindole with {methyl_count} methyl group(s)"
# Pr=0.675
# Recall=1.0