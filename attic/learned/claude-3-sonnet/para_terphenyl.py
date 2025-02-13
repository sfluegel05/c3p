"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl compounds
A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    Para-terphenyls have a central benzene ring with two phenyl groups attached in para positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for presence of at least 3 benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_matches) < 3:
        return False, "Less than 3 benzene rings found"
    
    # Look for para-terphenyl core structure
    # This pattern looks for a central benzene ring with two carbons in para positions
    # that are connected to other aromatic rings
    para_terphenyl_pattern = Chem.MolFromSmarts("c1([#6]c2ccccc2)cc([#6]c3ccccc3)ccc1")
    
    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        # Try alternative pattern that allows for more substitution
        alt_pattern = Chem.MolFromSmarts("c1(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)c([*])c([*])c(-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3)c([*])c1([*])")
        if not mol.HasSubstructMatch(alt_pattern):
            return False, "No para-terphenyl core structure found"
    
    # Additional validation - check ring count and connectivity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"
        
    # Count number of aromatic rings
    aromatic_rings = len([ring for ring in ring_info.AtomRings() 
                         if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)])
    if aromatic_rings < 3:
        return False, "Insufficient number of aromatic rings"

    return True, "Contains para-terphenyl core structure with appropriate substitution pattern"