"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:xxxxx lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    LPS typically contains a lipid A region (with phosphorylated glucosamine disaccharide and fatty acids),
    a core oligosaccharide, and an O-antigen polysaccharide.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for lipid A components: phosphorylated glucosamine disaccharide
    # Look for glucosamine (NH group in a sugar ring) and phosphate groups
    glucosamine_pattern = Chem.MolFromSmarts("[NH][C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)CO")
    if not mol.HasSubstructMatch(glucosamine_pattern):
        return False, "No glucosamine moiety (lipid A component)"
    
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups (lipid A component)"
    
    # Check for 3-hydroxytetradecanoic acid (14 carbons, hydroxyl at position 3)
    # Approximate pattern: 14 carbons in chain with OH on third carbon
    # SMARTS: [CH3]-[CH2]10-[CH](O)-[CH2]-[C](=O)-O-
    hydroxy_fa_pattern = Chem.MolFromSmarts("[CH2](-[CH2])-[CH](O)-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[C](=O)O")
    if not mol.HasSubstructMatch(hydroxy_fa_pattern):
        return False, "No 3-hydroxytetradecanoic acid units"
    
    # Check for oligosaccharide components (at least three sugar rings)
    sugar_ring_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)CO")
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if len(sugar_matches) < 3:
        return False, f"Only {len(sugar_matches)} sugar units, need at least 3"
    
    # Check for Kdo (octulosonic acid) component: carboxylic acid and ketone in a sugar-like structure
    kdo_pattern = Chem.MolFromSmarts("[C](=O)-C(=O)-O-[C@H]1O[C@@H]([C@@H](O)[C@H](O)[C@@H]1O)CO")
    if not mol.HasSubstructMatch(kdo_pattern):
        return False, "No Kdo (octulosonic acid) component"
    
    return True, "Contains lipid A, core oligosaccharide, and Kdo components"