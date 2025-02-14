"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:51139 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid carrying a hydroxy group on the
    carbon atom at position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for dicarboxylic acid pattern (two -C(=O)O groups)
    dicarboxylic_acid_pattern = Chem.MolFromSmarts("[$(C(=O)O)][$(C(=O)O)]")
    dicarboxylic_acid_matches = mol.GetSubstructMatches(dicarboxylic_acid_pattern)
    
    if not dicarboxylic_acid_matches:
        return False, "Not a dicarboxylic acid"
    
    # Find alpha carbon (between two carboxylic acids)
    alpha_carbon_idx = None
    for match in dicarboxylic_acid_matches:
        carboxyl_atom1 = mol.GetAtomWithIdx(match[0])
        carboxyl_atom2 = mol.GetAtomWithIdx(match[1])
        
        for bond in carboxyl_atom1.GetBonds():
            if bond.GetBeginAtomIdx() == carboxyl_atom2.GetIdx():
                continue
            alpha_carbon_idx = bond.GetBeginAtomIdx()
            break
        
        if alpha_carbon_idx is not None:
            break
    
    if alpha_carbon_idx is None:
        return False, "No alpha carbon found between two carboxylic acids"
    
    # Check for hydroxy group on alpha carbon
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    has_hydroxyl = any(neighbor.GetAtomicNum() == 8 and
                       neighbor.GetTotalNumHs() == 1
                       for neighbor in alpha_carbon.GetNeighbors())
    
    if not has_hydroxyl:
        return False, "No hydroxyl group on alpha carbon"
    
    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 300:
        return False, "Molecular weight outside typical range for 2-hydroxydicarboxylic acids"
    
    # Handle stereochemistry
    AllChem.AdjustQueryParameters()
    query_mol = Chem.MolFromSmarts("[C@@](C(=O)O)(O)[C@@](C(=O)O)(O)")
    if not mol.HasSubstructMatch(query_mol):
        return False, "Stereochemistry does not match 2-hydroxydicarboxylic acid pattern"
    
    return True, "Dicarboxylic acid with hydroxyl group on alpha carbon"