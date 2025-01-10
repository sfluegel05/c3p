"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: flavin compounds
Based on dimethylisoalloxazine skeleton with N10 substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavins are derivatives of dimethylisoalloxazine with N10 substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core tricyclic system with N10 position marked
    # [#7:1] marks the N10 position we want to check for substitution
    core_pattern = Chem.MolFromSmarts("""
        [#7:1]1-[#6]2=[#6][#6]=[#6][#6](=[#6]2)C
        :[#7][#6]3=[#7][#6](=[O])[#7][#6](=[O])[#6]3=[#7]1
    """)
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core tricyclic flavin system"

    # Get matches and check N10 substitution
    matches = mol.GetSubstructMatches(core_pattern)
    for match in matches:
        n10_idx = match[0]  # First atom in SMARTS pattern is our N10
        n10_atom = mol.GetAtomWithIdx(n10_idx)
        
        # Check degree of N10 (should be 3 for substitution)
        if n10_atom.GetDegree() < 3:
            return False, "N10 position not substituted"
            
    # Check for two methyl groups pattern
    dimethyl_pattern = Chem.MolFromSmarts("C-c1c(C)cc2")
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Missing characteristic methyl groups"

    # Verify presence of the two carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("[#6](=[O])[#7][#6](=[O])")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing characteristic carbonyl groups"

    return True, "Contains dimethylisoalloxazine core with N10 substitution"

def test_flavins():
    """Test function with known flavins"""
    test_cases = [
        # Riboflavin (positive case)
        "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1",
        # FMN (positive case)
        "C12=NC(NC(C1=NC=3C(N2C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O)=CC(=C(C3)C)C)=O)=O",
        # Lumiflavin (positive case)
        "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",
        # Not a flavin (benzene - negative case)
        "c1ccccc1",
        # Not a flavin (purine - negative case)
        "n1c2[nH]cnc2nc1"
    ]
    
    for smiles in test_cases:
        result, reason = is_flavin(smiles)
        print(f"SMILES: {smiles}")
        print(f"Is flavin: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_flavins()