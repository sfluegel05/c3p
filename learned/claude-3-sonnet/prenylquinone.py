"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
A quinone substituted by a polyprenyl-derived side-chain.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for quinone core (cyclohexa-2,5-diene-1,4-dione)
    quinone_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")  # para-quinone
    quinone_pattern2 = Chem.MolFromSmarts("C1(=O)C(=O)C=CC=C1") # ortho-quinone
    
    if not (mol.HasSubstructMatch(quinone_pattern) or mol.HasSubstructMatch(quinone_pattern2)):
        # Check for more complex quinone patterns that might have substituents
        complex_quinone = Chem.MolFromSmarts("C1(=O)C(~*)=C(~*)C(=O)C(~*)=C1~*")
        if not mol.HasSubstructMatch(complex_quinone):
            return False, "No quinone core structure found"

    # Look for prenyl/isoprenoid patterns
    prenyl_pattern = Chem.MolFromSmarts("CC(=C)C") # Basic prenyl group
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)CCC") # Isoprene unit
    
    if not (mol.HasSubstructMatch(prenyl_pattern) or mol.HasSubstructMatch(isoprene_pattern)):
        return False, "No prenyl/isoprenoid group found"
    
    # Check for common prenylquinone features
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:  # Minimum size for a basic prenylquinone
        return False, "Too few carbons for a prenylquinone"
        
    if o_count < 2:  # Need at least the two quinone oxygens
        return False, "Too few oxygens for a quinone structure"

    # Look for common substituent patterns
    methoxy_pattern = Chem.MolFromSmarts("COC")
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    
    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    if ring_count < 1:
        return False, "No rings found"
    
    # Look for long carbon chains (typical in prenylquinones)
    long_chain = Chem.MolFromSmarts("CCCCCC")
    has_long_chain = mol.HasSubstructMatch(long_chain)
    
    # Calculate number of rotatable bonds to help identify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Specific subclass patterns
    menaquinone_pattern = Chem.MolFromSmarts("C1=CC=CC2=C1C(=O)C(C)=C(CC=C)C2=O")
    ubiquinone_pattern = Chem.MolFromSmarts("COC1=C(OC)C(=O)C(C)=C(CC=C)C1=O")
    
    is_menaquinone = mol.HasSubstructMatch(menaquinone_pattern)
    is_ubiquinone = mol.HasSubstructMatch(ubiquinone_pattern)
    
    # Classify based on features
    features = []
    if has_methoxy:
        features.append("methoxy substituted")
    if has_long_chain:
        features.append("long chain")
    if is_menaquinone:
        features.append("menaquinone-like")
    if is_ubiquinone:
        features.append("ubiquinone-like")
    if n_rotatable >= 5:
        features.append("flexible side chain")
        
    if len(features) >= 2:  # If we have multiple prenylquinone-like features
        return True, f"Prenylquinone with features: {', '.join(features)}"
    elif n_rotatable >= 3 and (has_methoxy or has_long_chain):
        return True, "Basic prenylquinone structure"
    else:
        return False, "Missing characteristic prenylquinone features"