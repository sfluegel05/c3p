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

    # Look for carboxylate pattern with adjacent carbon
    # [O-]C(=O)C pattern where the second C can have other connections
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)C")
    
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"
    
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    for match in matches:
        carboxylate_o = match[0]  # [O-]
        carboxylate_c = match[1]  # C of C(=O)[O-]
        alpha_c = match[2]        # Adjacent C
        
        # Get the atoms
        carb_atom = mol.GetAtomWithIdx(carboxylate_c)
        alpha_atom = mol.GetAtomWithIdx(alpha_c)
        
        # Verify carboxylate carbon has exactly one [O-] and one =O
        o_minus_count = sum(1 for n in carb_atom.GetNeighbors() 
                          if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1)
        o_double_count = sum(1 for n in carb_atom.GetNeighbors() 
                           if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(carb_atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE)
        
        if o_minus_count != 1 or o_double_count != 1:
            continue
            
        # Look for oxo group on alpha carbon
        alpha_neighbors = alpha_atom.GetNeighbors()
        has_oxo = False
        for n in alpha_neighbors:
            if (n.GetAtomicNum() == 8 and 
                mol.GetBondBetweenAtoms(alpha_atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE):
                has_oxo = True
                break
                
        if not has_oxo:
            continue
            
        # Check that alpha carbon doesn't have another carboxylate group
        if sum(1 for n in alpha_atom.GetNeighbors() 
               if n.GetAtomicNum() == 6 and any(nn.GetFormalCharge() == -1 for nn in n.GetNeighbors())) > 0:
            continue
            
        # Count total number of carboxylate groups
        carboxylate_pattern_full = Chem.MolFromSmarts("[O-]C(=O)")
        if len(mol.GetSubstructMatches(carboxylate_pattern_full)) > 1:
            continue  # More than one carboxylate group
            
        # Check that this isn't part of an acid anhydride
        anhydride_pattern = Chem.MolFromSmarts("C(=O)OC(=O)")
        if mol.HasSubstructMatch(anhydride_pattern):
            continue
            
        return True, "Contains carboxylate with oxo group at 2-position"

    return False, "No valid 2-oxo carboxylate group found"