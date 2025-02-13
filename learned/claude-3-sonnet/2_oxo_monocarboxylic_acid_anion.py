"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:38478 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    Must have a carboxylate group with an oxo group at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate-oxo pattern
    # Match [O-]C(=O)CC pattern where C is the carboxylate carbon
    # and the second carbon has an oxo group
    pattern = Chem.MolFromSmarts("[O-]C(=O)C(=O)")
    
    if not mol.HasSubstructMatch(pattern):
        return False, "No 2-oxo carboxylate group found"
    
    matches = mol.GetSubstructMatches(pattern)
    
    for match in matches:
        carboxylate_o = match[0]  # [O-]
        carboxylate_c = match[1]  # C of C(=O)[O-]
        alpha_c = match[2]        # C of C(=O)R
        
        # Get the carboxylate carbon atom
        carb_atom = mol.GetAtomWithIdx(carboxylate_c)
        alpha_atom = mol.GetAtomWithIdx(alpha_c)
        
        # Verify carboxylate carbon has exactly one [O-] and one =O
        o_minus_count = sum(1 for n in carb_atom.GetNeighbors() 
                          if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1)
        o_double_count = sum(1 for n in carb_atom.GetNeighbors() 
                           if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(carb_atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE)
        
        if o_minus_count != 1 or o_double_count != 1:
            continue
            
        # Verify alpha carbon has exactly one =O
        alpha_o_double_count = sum(1 for n in alpha_atom.GetNeighbors() 
                                 if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(alpha_atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE)
        
        if alpha_o_double_count != 1:
            continue
            
        # Check that this isn't part of an acid anhydride
        anhydride_pattern = Chem.MolFromSmarts("C(=O)OC(=O)")
        if mol.HasSubstructMatch(anhydride_pattern):
            continue
            
        # Check that we don't have two carboxylate groups on the alpha carbon
        if sum(1 for n in alpha_atom.GetNeighbors() if n.GetAtomicNum() == 6 and any(nn.GetFormalCharge() == -1 for nn in n.GetNeighbors())) > 1:
            continue
            
        return True, "Contains carboxylate with oxo group at 2-position"

    return False, "No valid 2-oxo carboxylate group found"