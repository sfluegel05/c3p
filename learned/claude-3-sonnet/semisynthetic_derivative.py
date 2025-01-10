"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:35620 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.GraphDescriptors import BertzCT

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is likely a semisynthetic derivative based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Calculate basic properties
    num_atoms = mol.GetNumAtoms()
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    complexity = BertzCT(mol)
    num_stereocenters = len(Chem.FindMolChiralCenters(mol))
    
    # Too small or simple to be semisynthetic
    if num_atoms < 15:
        return False, "Too small to be a semisynthetic derivative"
    if num_rings < 2:
        return False, "Too few rings for a typical semisynthetic derivative"
    
    # Check for common modification patterns
    
    # Halogenation pattern
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    has_halogens = mol.HasSubstructMatch(halogen_pattern)
    
    # N-alkylation pattern
    n_alkyl_pattern = Chem.MolFromSmarts("[CH3,CH2][N;!$(N=*)]")
    has_n_alkyl = mol.HasSubstructMatch(n_alkyl_pattern)
    
    # O-acetylation pattern
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)O")
    has_acetyl = mol.HasSubstructMatch(acetyl_pattern)
    
    # O-methylation pattern
    methoxy_pattern = Chem.MolFromSmarts("CO")
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    
    # Phosphorylation pattern
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-,OH])([O-,OH])")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    
    # Count modifications
    num_modifications = sum([has_halogens, has_n_alkyl, has_acetyl, 
                           has_methoxy, has_phosphate])
    
    # Complex natural product-like features
    has_complex_ring_system = num_rings >= 3 and complexity > 500
    has_significant_chirality = num_stereocenters >= 2
    
    # Decision making
    if has_complex_ring_system and has_significant_chirality:
        if num_modifications >= 1:
            return True, "Complex natural product scaffold with synthetic modifications"
            
    if complexity > 1000 and num_modifications >= 2:
        return True, "Highly complex molecule with multiple synthetic modifications"
        
    if num_stereocenters >= 4 and num_modifications >= 1:
        return True, "Stereochemically complex molecule with synthetic modifications"
        
    if has_complex_ring_system and (has_halogens or has_phosphate):
        return True, "Natural product-like structure with characteristic synthetic modifications"
        
    if num_rings >= 4 and (has_n_alkyl or has_acetyl or has_methoxy):
        return True, "Multi-ring system with typical semisynthetic modifications"
        
    # Special cases for common semisynthetic classes
    penicillin_pattern = Chem.MolFromSmarts("S1CC2N(C1)C(=O)C2")
    if mol.HasSubstructMatch(penicillin_pattern):
        return True, "Modified β-lactam antibiotic structure"
        
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2(C1))")
    if mol.HasSubstructMatch(steroid_pattern) and num_modifications >= 1:
        return True, "Modified steroid structure"
        
    macrolide_pattern = Chem.MolFromSmarts("O=C1CCC(OC)CCC1")
    if mol.HasSubstructMatch(macrolide_pattern) and num_modifications >= 1:
        return True, "Modified macrolide structure"
    
    # If none of the above conditions are met
    if complexity < 300 or num_stereocenters < 2:
        return False, "Lacks complexity typical of semisynthetic derivatives"
        
    return False, "Does not match typical semisynthetic patterns"