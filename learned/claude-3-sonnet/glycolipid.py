"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_glycolipid, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sugar/carbohydrate patterns
    # Pattern for pyranose ring with multiple OH groups
    sugar_pattern = Chem.MolFromSmarts("[C]1[O][C]([C,O])[C]([O,N])[C]([O,N])[C]1[O,N]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    
    if not sugar_matches:
        return False, "No sugar moiety found"
    
    # Look for long carbon chains (fatty acids/lipid parts)
    fatty_chain_pattern = Chem.MolFromSmarts("[C][C][C][C][C][C][C][C]")
    fatty_chain_matches = mol.GetSubstructMatches(fatty_chain_pattern)
    
    if not fatty_chain_matches:
        return False, "No fatty acid/lipid chains found"
    
    # Look for glycosidic linkage (C-O-C connecting sugar)
    glycosidic_pattern = Chem.MolFromSmarts("[C][O][C]1[O][C]([C,O])[C][C][C]1[O,N]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
    
    # Additional checks for common glycolipid features
    
    # Count oxygens (should be multiple due to sugar hydroxyls)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Too few oxygen atoms for a glycolipid"
    
    # Check molecular weight (should be substantial due to sugar + lipid parts)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for typical glycolipid"
        
    # Look for common glycolipid features
    features = []
    
    # Check for sphingosine base
    sphingosine_pattern = Chem.MolFromSmarts("[CH2][CH]([NH])[CH]([OH])")
    if mol.HasSubstructMatch(sphingosine_pattern):
        features.append("sphingosine base")
    
    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if mol.HasSubstructMatch(glycerol_pattern):
        features.append("glycerol backbone")
    
    # Check for amide linkage (common in ceramide-based glycolipids)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if mol.HasSubstructMatch(amide_pattern):
        features.append("amide linkage")
    
    # Count sugar rings
    num_sugars = len(sugar_matches)
    sugar_desc = f"{num_sugars} sugar ring{'s' if num_sugars > 1 else ''}"
    features.append(sugar_desc)
    
    # Count fatty chains
    num_chains = len(fatty_chain_matches)
    chain_desc = f"{num_chains} fatty chain{'s' if num_chains > 1 else ''}"
    features.append(chain_desc)
    
    features_str = ", ".join(features)
    return True, f"Glycolipid containing {features_str}"