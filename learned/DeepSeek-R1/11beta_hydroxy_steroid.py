"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Basic steroid nucleus pattern (simplified four fused rings)
    steroid_smarts = Chem.MolFromSmarts("""
        C@H12CC[C@H]3C(C1)CC[C@H]4[C@@]3(C)CC[C@H]([C@@H]5CCCC(C)(C)C5)C4""")
    if not mol.HasSubstructMatch(steroid_smarts):
        return False, "No steroid nucleus found"
    
    # Pattern for 11beta-hydroxy group: hydroxyl at C11 with beta configuration
    # Beta configuration typically means same face as C18/C19 methyl groups
    hydroxy_smarts = Chem.MolFromSmarts("[C@@H](O)[C@H]1C[C@H]2CC[C@H]3C(C2)CC[C@H]4[C@@]3(C)CCC14")
    if not mol.HasSubstructMatch(hydroxy_smarts):
        return False, "No 11beta-hydroxy group detected"
    
    # Verify hydroxyl group exists at matched position
    matches = mol.GetSubstructMatches(hydroxy_smarts)
    for match in matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:
                return True, "11beta-hydroxy group on steroid nucleus"
    
    return False, "Missing hydroxyl group at position 11"