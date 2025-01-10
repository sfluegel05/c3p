"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Essential structural components of CoA
    
    # Adenine base (multiple patterns to catch different representations)
    adenine_patterns = [
        "[#7]1:[#6]:[#7]:[#6]2:[#6](:[#7]:[#6]:[#7]:[#6]:2:[#7]1)[#7]",  # Most generic form
        "n1cnc2c(N)ncnc12",  # Common form
        "N1C=NC2=C1N=CN=C2N"  # Alternative form
    ]
    
    has_adenine = False
    for pattern in adenine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_adenine = True
            break
    if not has_adenine:
        return False, "Missing adenine base structure"

    # Ribose with phosphate
    ribose_phosphate = Chem.MolFromSmarts("[CH2]OP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose phosphate group"

    # Pantetheine arm patterns (more flexible)
    pantetheine_patterns = [
        "NCCC(=O)NCCSC(=O)",  # Basic pattern
        "NC(=O)CCNC(=O)CCSC(=O)",  # Alternative pattern
        "[NH2][CH2][CH2][C](=[O])[NH][CH2][CH2]SC(=O)"  # More explicit pattern
    ]
    
    has_pantetheine = False
    for pattern in pantetheine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_pantetheine = True
            break
    if not has_pantetheine:
        return False, "Missing pantetheine arm structure"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester linkage"

    # Find the acyl chain
    acyl_carbon = thioester_matches[0][0]
    
    # Get the fragment containing the acyl chain
    # Break the S-C bond of thioester to isolate acyl part
    fragments = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(thioester_matches[0][0], 
                                         thioester_matches[0][1]).GetIdx()], addDummies=False)
    
    if fragments is None:
        return False, "Could not analyze acyl chain"
        
    # Count carbons in the acyl fragment
    acyl_fragment = Chem.GetMolFrags(fragments, asMols=True)[0]
    carbon_count = sum(1 for atom in acyl_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count > 6:
        return False, f"Acyl chain too long ({carbon_count} carbons) for short-chain fatty acid"
    if carbon_count < 2:
        return False, f"Acyl chain too short ({carbon_count} carbons)"

    # Check for modifications
    modifications = []
    mod_patterns = {
        "hydroxy": ("[OX2H1]", "hydroxyl"),
        "amino": ("[NX3H2]", "amino"),
        "unsaturated": ("[CX3]=[CX3]", "unsaturated"),
        "branched": ("[CX4]([CX4])([CX4])[CX4]", "branched")
    }
    
    for mod_name, (pattern, desc) in mod_patterns.items():
        if acyl_fragment.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            modifications.append(desc)

    mod_text = " with " + ", ".join(modifications) if modifications else ""
    return True, f"Short-chain ({carbon_count} carbons) fatty acyl-CoA{mod_text}"