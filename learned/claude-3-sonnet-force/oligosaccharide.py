"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:36973 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound where monosaccharide units are joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for monosaccharide rings using SMARTS pattern
    monosaccharide_pattern = Chem.MolFromSmarts("[OC1C(O)C(O)C(O)C(O)C1]")
    monosaccharide_rings = mol.GetSubstructMatches(monosaccharide_pattern)
    
    if not monosaccharide_rings:
        return False, "No monosaccharide rings found"
    
    # Look for glycosidic linkages (acetal bonds between rings) using SMARTS pattern
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[OC1C(O)C(O)C(O[C2C(O)C(O)C(O)C(O)C2])C(O)C1]")
    glycosidic_linkages = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    
    if not glycosidic_linkages:
        return False, "No glycosidic linkages found"
    
    # Count rotatable bonds (oligosaccharides typically have few rotatable bonds)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_monosaccharide_rings = len(monosaccharide_rings)
    n_glycosidic_linkages = len(glycosidic_linkages)
    
    # Adjust rotatable bond threshold based on the size of the oligosaccharide
    rotatable_bond_threshold = 10 + (n_monosaccharide_rings + n_glycosidic_linkages - 2) * 2
    
    if n_rotatable > rotatable_bond_threshold:
        return False, f"Too many rotatable bonds ({n_rotatable}) for an oligosaccharide with {n_monosaccharide_rings} monosaccharide rings and {n_glycosidic_linkages} glycosidic linkages"
    
    return True, f"Contains {n_monosaccharide_rings} monosaccharide rings connected by {n_glycosidic_linkages} glycosidic linkages"