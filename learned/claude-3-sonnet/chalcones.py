"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

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
    
    # Basic chalcone patterns with more flexible matching
    # Core pattern including variations in aromatic systems
    chalcone_patterns = [
        # Classic chalcone pattern
        Chem.MolFromSmarts("[$(a:a)]!@[CH]=[CH]!@C(=O)!@[$(a:a)]"),
        # Dihydrochalcone pattern
        Chem.MolFromSmarts("[$(a:a)]!@[CH2][CH2]!@C(=O)!@[$(a:a)]"),
        # Pattern for fused ring systems
        Chem.MolFromSmarts("[$(a:a:a)]!@[CH]=[CH]!@C(=O)!@[$(a:a)]"),
        # Pattern allowing for substituted double bond
        Chem.MolFromSmarts("[$(a:a)]!@[C]=[C]!@C(=O)!@[$(a:a)]")
    ]
    
    found_pattern = False
    for pattern in chalcone_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_pattern = True
            break
            
    if not found_pattern:
        return False, "No chalcone core structure found"

    # Check for common substitution patterns
    substitution_patterns = {
        "hydroxy": Chem.MolFromSmarts("cO"),
        "methoxy": Chem.MolFromSmarts("cOC"),
        "prenyl": Chem.MolFromSmarts("CC(C)=CC"),
        "methylenedioxy": Chem.MolFromSmarts("OCO"),
    }
    
    # Count substitutions
    substitutions = []
    for name, pattern in substitution_patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            substitutions.append(name)
    
    # Count key features
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    ketone_count = len(mol.GetSubstructMatches(ketone_pattern)) if ketone_pattern else 0
    
    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 2:
        return False, "Insufficient number of aromatic rings"

    # Check molecular properties
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 900:  # Adjusted range based on examples
        return False, "Molecular weight outside typical range for chalcones"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count == 0:
        return False, "No oxygen atoms found"
    
    if c_count < 15:  # Minimum carbons for basic chalcone structure
        return False, "Too few carbons for chalcone structure"

    # Build classification reason
    substitution_desc = ""
    if substitutions:
        substitution_desc = f" with {', '.join(substitutions)} groups"
    
    chalcone_type = "dihydrochalcone" if mol.HasSubstructMatch(chalcone_patterns[1]) else "chalcone"
    reason = f"Contains {chalcone_type} core structure{substitution_desc}"

    return True, reason