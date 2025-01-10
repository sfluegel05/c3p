"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a 1,3-diphenylpropenone (benzylideneacetophenone) or its derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic chalcone pattern: Ar-CH=CH-C(=O)-Ar
    # Looking for:
    # - Two aromatic rings
    # - Connected by a propenone bridge (-CH=CH-C(=O)-)
    chalcone_pattern = Chem.MolFromSmarts("[a;r6]!@[CH]=[CH]!@C(=O)!@[a;r6]")
    
    # Alternative pattern for dihydrochalcones (reduced double bond)
    dihydrochalcone_pattern = Chem.MolFromSmarts("[a;r6]!@[CH2][CH2]!@C(=O)!@[a;r6]")
    
    matches = mol.GetSubstructMatches(chalcone_pattern)
    dihydro_matches = mol.GetSubstructMatches(dihydrochalcone_pattern)
    
    if not (matches or dihydro_matches):
        return False, "No chalcone or dihydrochalcone core structure found"

    # Count number of aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 2:
        return False, "Insufficient number of aromatic rings"

    # Additional checks for common chalcone features
    # Check for common substituents (OH, OMe, etc.)
    substituents = {
        "hydroxy": Chem.MolFromSmarts("[OH]"),
        "methoxy": Chem.MolFromSmarts("[OH0][CH3]"),
    }
    
    has_substituents = any(
        len(mol.GetSubstructMatches(pattern)) > 0 
        for pattern in substituents.values()
    )

    # Determine if it's a regular chalcone or dihydrochalcone
    if matches:
        chalcone_type = "chalcone"
    else:
        chalcone_type = "dihydrochalcone"

    # Build detailed classification reason
    substituent_desc = " with common substituents" if has_substituents else ""
    reason = f"Contains {chalcone_type} core structure (1,3-diphenylpropenone){substituent_desc}"
    
    # Additional validation to reduce false positives
    # Check molecular weight (most chalcones are between 200-600 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for chalcones"

    # Check carbon to oxygen ratio (typically between 3:1 and 10:1)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count == 0 or c_count/o_count > 10 or c_count/o_count < 3:
        return False, "Carbon to oxygen ratio outside typical range for chalcones"

    return True, reason