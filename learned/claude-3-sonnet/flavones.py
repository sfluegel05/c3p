"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: flavones
A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic flavone core structure (2-phenylchromen-4-one)
    # Using a more flexible SMARTS pattern that accounts for aromaticity
    # and allows for substitutions
    flavone_core = Chem.MolFromSmarts("c1c2c(oc(-c3ccccc3)cc2=O)cc1")
    
    if not mol.HasSubstructMatch(flavone_core):
        # Try alternative SMARTS pattern that's even more flexible
        alt_flavone_core = Chem.MolFromSmarts("c1c2c(oc(-[c;r6][c,n][c,n][c,n][c,n][c,n]3)cc2=O)cc1")
        if not mol.HasSubstructMatch(alt_flavone_core):
            return False, "Missing flavone core structure (2-phenylchromen-4-one skeleton)"
    
    # Verify the presence of key structural features
    
    # Check for the chromone (benzopyran-4-one) system
    chromone = Chem.MolFromSmarts("c1c2c(occ2=O)cc1")
    if not mol.HasSubstructMatch(chromone):
        return False, "Missing chromone (benzopyran-4-one) system"
    
    # Count rings to ensure we have at least two aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    if aromatic_rings < 2:
        return False, "Must contain at least two aromatic rings"
    
    # Check for minimum number of carbons and oxygens
    # Basic flavone should have at least 15 carbons and 2 oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for flavone structure"
    if o_count < 2:
        return False, f"Too few oxygens ({o_count}) for flavone structure"
    
    # Look for common substitution patterns
    substitutions = []
    
    # Hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("cO")
    if mol.HasSubstructMatch(hydroxy_pattern):
        substitutions.append("hydroxy")
    
    # Methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("cOC")
    if mol.HasSubstructMatch(methoxy_pattern):
        substitutions.append("methoxy")
    
    # Glycoside patterns
    glycoside_pattern = Chem.MolFromSmarts("[OH1,OH0][CH1][OH1,OH0][CH1][OH1,OH0][CH1,CH2]O")
    if mol.HasSubstructMatch(glycoside_pattern):
        substitutions.append("glycoside")
    
    # Prenyl groups
    prenyl_pattern = Chem.MolFromSmarts("CC(C)=CC")
    if mol.HasSubstructMatch(prenyl_pattern):
        substitutions.append("prenyl")
        
    substitution_str = ""
    if substitutions:
        substitution_str = f" with {', '.join(substitutions)} substitutions"
    
    return True, f"Contains flavone core structure (2-phenylchromen-4-one){substitution_str}"