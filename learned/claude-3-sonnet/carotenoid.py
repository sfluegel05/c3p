"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23042 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons - carotenoids typically have ~40 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:  # Allow some flexibility for degraded carotenoids
        return False, f"Too few carbons ({c_count}) for a carotenoid"
    
    # Look for polyene chain pattern (alternating single-double bonds)
    polyene_pattern = Chem.MolFromSmarts("C=CC=CC=CC=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No characteristic polyene chain found"
    
    # Count conjugated double bonds
    conjugated_matches = len(mol.GetSubstructMatches(polyene_pattern))
    if conjugated_matches < 2:
        return False, "Insufficient conjugated double bond system"
    
    # Look for common end groups
    cyclohexene_pattern = Chem.MolFromSmarts("C1C=C(C)CCC1")  # Beta-type end group
    cyclopentene_pattern = Chem.MolFromSmarts("C1C=C(C)CC1")   # Epsilon-type end group
    acyclic_pattern = Chem.MolFromSmarts("CC(C)=CCC")         # Acyclic end group
    
    has_typical_end = (mol.HasSubstructMatch(cyclohexene_pattern) or 
                      mol.HasSubstructMatch(cyclopentene_pattern) or
                      mol.HasSubstructMatch(acyclic_pattern))
    
    if not has_typical_end:
        return False, "No characteristic end groups found"
    
    # Check for common modifications
    hydroxy_pattern = Chem.MolFromSmarts("CO")  # xanthophylls
    keto_pattern = Chem.MolFromSmarts("CC(=O)C")
    epoxy_pattern = Chem.MolFromSmarts("COC")
    
    modifications = []
    if mol.HasSubstructMatch(hydroxy_pattern):
        modifications.append("hydroxylated")
    if mol.HasSubstructMatch(keto_pattern):
        modifications.append("keto")
    if mol.HasSubstructMatch(epoxy_pattern):
        modifications.append("epoxidated")
        
    # Calculate degree of unsaturation
    num_double_bonds = rdMolDescriptors.CalcNumBonds(mol) - mol.GetNumAtoms() + 1
    
    # Final checks and classification
    if num_double_bonds < 8:
        return False, "Insufficient unsaturation for a carotenoid"
        
    mod_str = " (" + ", ".join(modifications) + ")" if modifications else ""
    return True, f"C{c_count} terpenoid with conjugated polyene system{mod_str}"