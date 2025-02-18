"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:79126 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative with oxygen/hydroxyl groups replaced by sulfur.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for sulfur presence
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Look for carbohydrate-like structure (pyranose/furanose with multiple hydroxyls)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O1)O"),  # pyranose template
        Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@@H](O1)O")           # furanose template
    ]
    sugar_matches = any(mol.HasSubstructMatch(patt) for patt in sugar_patterns if patt)
    
    # If no classic sugar pattern, check for modified structures with >=2 hydroxyls
    if not sugar_matches:
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() 
                            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
        if hydroxyl_count < 2:
            return False, "Insufficient hydroxyl groups for carbohydrate derivative"

    # Check sulfur positions relative to carbohydrate structure
    for sulfur in sulfur_atoms:
        # Case 1: Sulfur in ring position (replacing oxygen)
        if sulfur.IsInRing():
            ring = next((r for r in mol.GetRingInfo().AtomRings() if sulfur.GetIdx() in r), None)
            if ring and len(ring) in (5,6):
                # Check if this could be a sugar ring with S replacing O
                return True, "Sulfur in carbohydrate ring position"
        
        # Case 2: Sulfur attached to ring carbon (replacing hydroxyl/oxygen)
        for neighbor in sulfur.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.IsInRing():
                ring = next((r for r in mol.GetRingInfo().AtomRings() if neighbor.GetIdx() in r), None)
                if ring and len(ring) in (5,6):
                    # Check if adjacent to oxygen-rich environment (carbohydrate-like)
                    o_count = sum(1 for n in neighbor.GetNeighbors() 
                                 if n.GetAtomicNum() == 8 and n.GetTotalNumHs() >= 1)
                    if o_count >= 1:
                        return True, "Sulfur attached to carbohydrate ring carbon"
        
        # Case 3: Sulfur in glycosidic linkage (thioether between two sugars)
        glycosidic_s = Chem.MolFromSmarts("[C][SX2][C@H]1O[C@H]")
        if mol.HasSubstructMatch(glycosidic_s):
            return True, "Sulfur in glycosidic linkage"

    return False, "No sulfur substitution in carbohydrate structure"