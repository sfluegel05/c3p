"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight (CoA itself is quite large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 750:  # CoA MW is ~767
        return False, "Molecular weight too low for acyl-CoA"

    # Count key atoms
    num_P = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15])  # Phosphorus
    num_S = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16])  # Sulfur
    num_N = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])   # Nitrogen
    
    # Basic CoA checks
    if num_P < 3:
        return False, "Missing phosphate groups characteristic of CoA"
    if num_S < 1:
        return False, "Missing sulfur atom required for thioester linkage"
    if num_N < 5:
        return False, "Insufficient nitrogen atoms for CoA structure"

    # Check for key CoA structural features
    coenzyme_a_patterns = [
        # Thioester linkage
        "[CX3](=O)[SX2]",
        # Adenine base
        "n1c(nc2c1ncnc2N)",
        # Pantetheine part
        "CC(C)(COP(=O)(O)OP(=O)(O)O)C"
    ]
    
    for pattern in coenzyme_a_patterns:
        if not mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Missing essential CoA structural feature"

    # Check for branching in the fatty acid portion
    branching_patterns = [
        # General branched carbon
        "[CH1,CH0](-[CH3])(-[CH2,CH3])-[!R]",
        # Terminal isopropyl group
        "[CH3]-[CH](-[CH3])-[CH2]",
        # Internal branching
        "[CH2]-[CH](-[CH3])-[CH2]"
    ]
    
    found_branching = False
    for pattern in branching_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_branching = True
            break
            
    if not found_branching:
        return False, "No branching pattern found in fatty acid portion"

    # Check for fatty acid chain length (at least 4 carbons)
    carbon_chain = "[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(carbon_chain)):
        return False, "Carbon chain too short for fatty acid"

    # Look for common modifications
    modifications = []
    
    # Check for hydroxy groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH1,CH0]([OH1])")):
        modifications.append("hydroxy-modified")
    
    # Check for unsaturation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        modifications.append("unsaturated")
    
    # Check for oxo groups (besides thioester)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]C(=O)[CH2]")):
        modifications.append("oxo-modified")

    base_reason = "Contains CoA moiety with thioester linkage and branched-chain fatty acid pattern"
    if modifications:
        return True, f"{base_reason} ({', '.join(modifications)})"
    return True, base_reason