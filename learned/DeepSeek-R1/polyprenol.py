"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:26188 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a prenol with the general formula H-[CH2C(Me)=CHCH2]nOH (n > 1),
    composed of multiple isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]-[C]'))
    if not hydroxyl_matches:
        return False, "No hydroxyl group present"
    
    # Verify hydroxyl is on a terminal carbon (carbon with degree 1)
    terminal_hydroxyl = False
    for match in hydroxyl_matches:
        carbon_idx = match[1]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        if carbon.GetDegree() == 1:
            terminal_hydroxyl = True
            break
    if not terminal_hydroxyl:
        return False, "Hydroxyl group not on terminal carbon"
    
    # Find methyl groups attached to double-bonded carbons (isoprene unit markers)
    # SMARTS: Carbon with a double bond and a methyl group
    isoprene_unit_pattern = Chem.MolFromSmarts('[CH0]([CH3])=*')
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit_pattern)
    
    # Need at least two isoprene units (n > 1)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units (requires >1)"
    
    return True, "Terminal hydroxyl with ≥2 isoprene units (methyl on double bond carbons)"