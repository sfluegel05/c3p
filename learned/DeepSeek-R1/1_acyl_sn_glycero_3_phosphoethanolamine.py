"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine (CHEBI:17517)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    Must have:
    - Glycerol backbone with (R)-configuration at sn-2
    - 1 acyl chain at sn-1 (ester linkage)
    - Phosphoethanolamine group at sn-3
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphoethanolamine group (O-P-O-C-C-N, allowing for charges)
    pepatt = Chem.MolFromSmarts("[PX4]=O.OCC[NH2,NH3+]")
    if not mol.HasSubstructMatch(pepatt):
        return False, "Missing phosphoethanolamine group"

    # Find glycerol backbone with correct substituents and stereochemistry
    # Pattern: sn-1: [CH2]OC(=O), sn-2: [C@H](O), sn-3: [CH2]OP(=O)...OCCN
    # Allow variations in phosphate protonation and ethanolamine charge
    core_patterns = [
        Chem.MolFromSmarts("[CH2]OC(=O)-[C@H](O)-[CH2]OP(=O)([O-])OCC[NH2,NH3+]"),
        Chem.MolFromSmarts("[CH2]OC(=O)-[C@H](O)-[CH2]OP(=O)(O)OCC[NH2,NH3+]"),
        # Alternative arrangement where phosphate is connected via oxygen first
        Chem.MolFromSmarts("O-P(=O)([O-])OCC[NH2,NH3+]-O-C-[C@H](O)-COC(=O)"),
    ]
    
    found_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            found_core = True
            break
    if not found_core:
        return False, "Missing core structure with R configuration"

    # Verify only one ester group (the 1-acyl)
    ester_patt = Chem.MolFromSmarts("[OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_patt)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester group, found {len(ester_matches)}"

    # Check acyl chain length (at least 10 carbons in the chain)
    # Find the ester oxygen, then trace the acyl chain
    ester_match = ester_matches[0]
    oxygen_idx = ester_match[0]
    oxygen = mol.GetAtomWithIdx(oxygen_idx)
    # The ester oxygen is connected to the carbonyl carbon (index 1 in ester_match)
    carbonyl_carbon = ester_match[1]
    # Traverse the chain starting from the carbonyl carbon
    chain_carbons = set()
    stack = [(mol.GetAtomWithIdx(carbonyl_carbon), 0)]  # (atom, depth)
    visited = set()
    while stack:
        atom, depth = stack.pop()
        if atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        if atom.GetAtomicNum() == 6 and atom.GetIdx() != carbonyl_carbon:
            chain_carbons.add(atom.GetIdx())
            # Continue traversing only carbons not part of the glycerol backbone
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in [carbonyl_carbon, oxygen.GetIdx()]:
                    stack.append((neighbor, depth + 1))
    if len(chain_carbons) < 10:
        return False, f"Acyl chain too short (only {len(chain_carbons)} carbons)"

    return True, "Has 1-acyl-sn-glycero-3-phosphoethanolamine structure"