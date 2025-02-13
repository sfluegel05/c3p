"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:36345 lipopolysaccharide 
Lipopolysaccharides are natural compounds consisting of a trisaccharide repeating unit 
(two heptose units and octulosonic acid) with oligosaccharide side chains and 
3-hydroxytetradecanoic acid units. They are a major constituent of the cell walls of Gram-negative bacteria.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for patterns of trisaccharide repeating unit
    trisaccharide_patterns = [
        Chem.MolFromSmarts("[C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)O[C@@H]2[C@@H]([C@H]([C@@H](O2)CO)O)O"),  # heptose-heptose-octulosonic acid
        Chem.MolFromSmarts("[C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)O[C@@H]2[C@@H]([C@H]([C@@H](O2)CO)O)O[C@H]3[C@@H]([C@H]([C@@H](O3)CO)O)O")  # heptose-heptose-octulosonic acid-heptose
    ]
    trisaccharide_match = any(mol.HasSubstructMatch(pattern) for pattern in trisaccharide_patterns)
    if not trisaccharide_match:
        return False, "No trisaccharide repeating unit found"
    
    # Look for oligosaccharide side chains
    oligosaccharide_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]~[OX2][CX4][OX2]~[OX2][CX4][OX2]")
    oligosaccharide_matches = mol.GetSubstructMatches(oligosaccharide_pattern)
    if not oligosaccharide_matches:
        return False, "No oligosaccharide side chains found"
    
    # Look for 3-hydroxytetradecanoic acid unit
    lipid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCC[C@@H](O)C(O)=O"),  # 3-hydroxytetradecanoic acid
        Chem.MolFromSmarts("[C@@H](CCCCCCCCCCCCC)O"),  # 3-hydroxylated lipid chain
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(O)=O")  # tetradecanoic acid
    ]
    lipid_match = any(mol.HasSubstructMatch(pattern) for pattern in lipid_patterns)
    if not lipid_match:
        return False, "No 3-hydroxytetradecanoic acid unit found"
    
    # Check molecular weight - lipopolysaccharides typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for lipopolysaccharide"
    
    return True, "Contains trisaccharide repeating unit, oligosaccharide side chains, and 3-hydroxytetradecanoic acid unit"