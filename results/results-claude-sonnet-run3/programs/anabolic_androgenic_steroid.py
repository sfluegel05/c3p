from rdkit import Chem
from rdkit.Chem import AllChem

def is_anabolic_androgenic_steroid(smiles: str):
    """
    Determines if a molecule is an anabolic androgenic steroid based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anabolic androgenic steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for steroid core structure using a more flexible SMARTS pattern
    # This matches the 4-ring system with more flexibility for different bond types
    steroid_core = '[#6]~1~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~1~[#6]~2~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~12'
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(steroid_core))
    if not matches:
        return False, "Does not contain steroid core structure"

    # Check for ketone or hydroxyl groups at typical positions
    c3_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]1~[#6]~[#6](=[O])~[#6]~[#6]~[#6]1'))
    c3_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]1~[#6]~[#6](O)~[#6]~[#6]~[#6]1'))
    c17_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]~[#6](=O)~[#6]'))
    c17_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]~[#6](O)~[#6]'))

    # Count carbon atoms to verify it's in the right range for steroids
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 18 or num_carbons > 25:
        return False, "Carbon count outside typical range for anabolic steroids"

    # Count rings
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 4:
        return False, "Does not contain minimum 4 rings required for steroid scaffold"

    # If molecule has the core structure and at least one ketone or hydroxyl group
    if (c3_ketone or c3_hydroxyl or c17_ketone or c17_hydroxyl):
        modifications = []
        if c3_ketone:
            modifications.append("C3 ketone")
        if c3_hydroxyl:
            modifications.append("C3 hydroxyl")
        if c17_ketone:
            modifications.append("C17 ketone")
        if c17_hydroxyl:
            modifications.append("C17 hydroxyl")
        
        return True, f"Anabolic androgenic steroid with {', '.join(modifications)}"
    
    return False, "Missing characteristic ketone or hydroxyl groups"
# Pr=None
# Recall=0.0