"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: isothiocyanate compounds (R-N=C=S)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates have the general formula R-N=C=S where R is an organic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N=C=S group
    isothiocyanate_pattern = Chem.MolFromSmarts("[NX2]=[CX2]=[SX1]")
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    
    if not matches:
        return False, "No isothiocyanate group (N=C=S) found"
    
    # For each potential match, verify:
    # 1. The nitrogen should have exactly one additional single bond (to R group)
    # 2. The carbon should only have the double bonds to N and S
    # 3. The sulfur should only have the double bond to C
    
    valid_groups = 0
    for match in matches:
        n_atom, c_atom, s_atom = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Check nitrogen connectivity
        if n_atom.GetDegree() != 2:  # One single bond to R, one double bond to C
            continue
            
        # Check carbon connectivity
        if c_atom.GetDegree() != 2:  # Should only connect to N and S
            continue
            
        # Check sulfur connectivity
        if s_atom.GetDegree() != 1:  # Should only connect to C
            continue
            
        valid_groups += 1
    
    if valid_groups == 0:
        return False, "No valid isothiocyanate groups found"
    
    # Check if there's at least one carbon or aromatic group attached
    # (to verify it's an organic compound)
    has_organic = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 or atom.GetIsAromatic():
            has_organic = True
            break
    
    if not has_organic:
        return False, "No organic group (R) found"
        
    if valid_groups == 1:
        return True, "Contains one isothiocyanate group (R-N=C=S)"
    else:
        return True, f"Contains {valid_groups} isothiocyanate groups"