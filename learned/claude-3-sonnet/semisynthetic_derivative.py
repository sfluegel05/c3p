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
    
    # Too small to be semisynthetic
    if num_atoms < 10:
        return False, "Too small to be a semisynthetic derivative"

    # Natural product core scaffolds (expanded)
    natural_scaffolds = {
        "steroid": "C1CC2CCC3C(C2(C1))",
        "beta_lactam": "N1C(=O)C2CSC21",
        "macrolide": "[#6]1[#6][#6][#6]([#8])[#6][#6][#6][#6]1",
        "tetracycline": "C1C(=O)C=C2C(=O)C3(O)C(=O)C(C(N)=O)C(=O)C3CC2C1",
        "alkaloid_1": "C1CN2CCC1CC2",
        "alkaloid_2": "c1ccc2c(c1)c1[nH]cc(C)c1n2", # β-carboline core
        "alkaloid_3": "C1CN2C(CC1)Cc1ccccc12", # ergoline core
        "flavonoid": "C1=CC(=O)C2=C(C=C(O)C=C2)O1",
        "morphinan": "C1CC2C3CC4(C2)C(C1)N(CC3)Cc1ccccc14",
        "cephalosporin": "S1CC2=C(N3C(=O)[C@H](NC)C3SC1)C(=O)O2",
        "penicillin": "S1CC2N(C(=O)C2)C(C(=O)O)C1"
    }

    # Synthetic modification patterns (refined)
    synthetic_mods = {
        "alkylation": "[CH3,CH2][N,O,S;!$(N=*);!$(S=O);!$(C(=O)O)]",
        "acylation": "[NX3,OX2][CX3](=O)[#6;!$(C=O)]",
        "halogenation": "[F,Cl,Br,I]",
        "phosphorylation": "[PX4](=O)([O-,OH])([O-,OH])",
        "sulfonation": "S(=O)(=O)[O-,OH]",
        "methylation": "[$([CH3][O,N,S]);!$(CC=O)]",
        "esterification": "[#6]C(=O)O[#6]",
        "amidation": "[NX3;H1,H0;$(NC=O)]"
    }

    # Check for natural product scaffolds
    scaffolds_found = []
    for name, pattern in natural_scaffolds.items():
        if check_smarts_match(mol, pattern):
            scaffolds_found.append(name)

    # Check for synthetic modifications
    mods_found = []
    for name, pattern in synthetic_mods.items():
        if check_smarts_match(mol, pattern):
            mods_found.append(name)

    # Carbohydrate pattern to exclude natural oligosaccharides
    carbohydrate_pattern = "[OX2H][CH]1[OX2][CH]([CH]([OX2H])[CH]([OX2H])[CH]1[OX2H])"
    is_carbohydrate = check_smarts_match(mol, carbohydrate_pattern)
    
    # Decision making
    if len(scaffolds_found) > 0:
        if len(mods_found) >= 1:
            return True, f"Natural product scaffold ({', '.join(scaffolds_found)}) with synthetic modifications ({', '.join(mods_found)})"
    
    # Handle special cases for alkaloid derivatives
    if "alkaloid_2" in scaffolds_found or "alkaloid_3" in scaffolds_found:
        if "alkylation" in mods_found or "acylation" in mods_found:
            return True, f"Modified alkaloid with {', '.join(mods_found)}"
    
    # Complex molecule with clear synthetic modifications
    if complexity > 600 and len(mods_found) >= 2:
        if not is_carbohydrate:
            return True, f"Complex molecule with synthetic modifications ({', '.join(mods_found)})"
    
    # Stereocomplex molecules with specific modifications
    if num_stereocenters >= 3 and ("halogenation" in mods_found or "phosphorylation" in mods_found):
        if not is_carbohydrate:
            return True, f"Stereocomplex molecule with specific modifications ({', '.join(mods_found)})"
    
    # Reject cases
    if is_carbohydrate and len(scaffolds_found) == 0:
        return False, "Natural oligosaccharide pattern"
    
    if complexity < 200 and len(scaffolds_found) == 0:
        return False, "Too simple for semisynthetic derivative"
    
    if num_rings == 0 and len(mods_found) < 2:
        return False, "Lacks characteristic ring systems and modifications"
    
    # Default case
    return False, "Does not match semisynthetic patterns"