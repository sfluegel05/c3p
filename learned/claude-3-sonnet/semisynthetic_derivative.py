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

def check_smarts_match(mol, smarts_pattern):
    """Helper function to safely check SMARTS matches"""
    pattern = Chem.MolFromSmarts(smarts_pattern)
    if pattern is None:
        return False
    return mol.HasSubstructMatch(pattern)

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
    if num_atoms < 12:
        return False, "Too small to be a semisynthetic derivative"
    
    # Check for common modification patterns
    modification_patterns = {
        "halogen": "[F,Cl,Br,I]",
        "n_alkyl": "[CH3,CH2][N;!$(N=*)]",
        "o_alkyl": "[CH3,CH2]O",
        "acetyl": "CC(=O)O",
        "methoxy": "COC",
        "phosphate": "[PX4](=O)([O-,OH])([O-,OH])",
        "sulfate": "OS(=O)(=O)[O-,OH]",
        "amine": "[NX3;H2,H1;!$(NC=O)]",
        "amide": "[NX3;H1,H0;$(NC=O)]",
        "ester": "[#6]C(=O)O[#6]"
    }
    
    modifications_found = []
    for mod_name, pattern in modification_patterns.items():
        if check_smarts_match(mol, pattern):
            modifications_found.append(mod_name)
    
    # Natural product scaffold patterns
    scaffold_patterns = {
        "steroid": "C1CC2CCC3C(C2(C1))",
        "beta_lactam": "N1C(=O)C2CSC21",
        "macrolide": "O=C1CCC(OC)CCC1",
        "tetracycline": "C1C(=O)C=C2C(=O)C3(O)C(=O)C(C(N)=O)C(=O)C3CC2C1",
        "alkaloid": "C1CN2CCC1CC2",
        "flavonoid": "C1=CC(=O)C2=C(C=C(O)C=C2)O1"
    }
    
    scaffolds_found = []
    for scaffold_name, pattern in scaffold_patterns.items():
        if check_smarts_match(mol, pattern):
            scaffolds_found.append(scaffold_name)
    
    # Decision making with adjusted criteria
    num_modifications = len(modifications_found)
    
    # Complex natural product with modifications
    if len(scaffolds_found) > 0 and num_modifications >= 1:
        return True, f"Natural product scaffold ({', '.join(scaffolds_found)}) with synthetic modifications ({', '.join(modifications_found)})"
    
    # High complexity with multiple modifications
    if complexity > 800 and num_modifications >= 2:
        return True, f"Complex molecule (complexity={int(complexity)}) with multiple modifications ({', '.join(modifications_found)})"
    
    # Significant stereochemistry with modifications
    if num_stereocenters >= 3 and num_modifications >= 1:
        return True, f"Stereochemically complex molecule ({num_stereocenters} stereocenters) with modifications ({', '.join(modifications_found)})"
    
    # Ring system with modifications
    if num_rings >= 3 and num_modifications >= 2:
        return True, f"Multi-ring system ({num_rings} rings) with synthetic modifications ({', '.join(modifications_found)})"
    
    # Special cases for specific modifications
    if "halogen" in modifications_found and complexity > 400:
        return True, "Halogenated complex molecule"
    
    if "phosphate" in modifications_found and num_rings >= 2:
        return True, "Phosphorylated ring system"
    
    if num_modifications >= 3 and num_atoms > 20:
        return True, f"Multiple synthetic modifications ({', '.join(modifications_found)})"
    
    # Reject cases
    if complexity < 250 or num_stereocenters == 0:
        return False, "Insufficient complexity for semisynthetic derivative"
    
    if num_rings == 0 and num_modifications < 2:
        return False, "Lacks characteristic ring systems and modifications"
    
    # Default case
    return False, "Does not match typical semisynthetic patterns"