"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid (CHEBI:24404)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid has a carbohydrate attached via glycosidic linkage to a sphingoid or ceramide.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for ceramide (amide group) or sphingoid (NH2 and adjacent OH)
    ceramide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[CX4]")
    sphingoid_pattern = Chem.MolFromSmarts("[NH2]-[CX4]-[CX4]-[OH]")
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)
    has_sphingoid = mol.HasSubstructMatch(sphingoid_pattern)
    
    if not (has_ceramide or has_sphingoid):
        return False, "No ceramide or sphingoid base found"
    
    # Check for glycosidic oxygen connected to a sugar (ether linkage to a ring with multiple OHs)
    # General sugar pattern: any ring with at least two hydroxyls
    sugar_pattern = Chem.MolFromSmarts("[C]1@[C]@[C]@[C]@[C]@[C]1(-[OH])-[OH]")
    glycosidic_o_pattern = Chem.MolFromSmarts("[OX2][C]@[C]")
    
    # Check if there's a glycosidic oxygen connected to a sugar
    has_glycosidic = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2:  # ether oxygen
            neighbors = atom.GetNeighbors()
            # Check if one side is part of ceramide/sphingoid and the other is a sugar
            # Simplification: check if sugar pattern exists and is connected via this oxygen
            # This is a heuristic and may not cover all cases
            for neighbor in neighbors:
                if neighbor.IsInRing() and mol.HasSubstructMatch(sugar_pattern):
                    has_glycosidic = True
                    break
            if has_glycosidic:
                break
    
    if not has_glycosidic:
        return False, "No glycosidic linkage to a sugar found"
    
    return True, "Contains carbohydrate attached to ceramide/sphingoid via glycosidic bond"