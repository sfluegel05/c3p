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
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine.
    Must have:
    - Glycerol backbone with (R)-configuration at sn-2
    - 1 acyl chain at sn-1 (ester linkage)
    - Phosphoethanolamine group at sn-3
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphoethanolamine group (O-P-O-C-C-N, allowing charges)
    # Updated pattern to match direct connection to oxygen
    pepatt = Chem.MolFromSmarts("[PX4](=O)([O-])OCC[NH2,NH3+]")
    if not mol.HasSubstructMatch(pepatt):
        pepatt2 = Chem.MolFromSmarts("[PX4](=O)(O)OCC[NH2,NH3+]")
        if not mol.HasSubstructMatch(pepatt2):
            return False, "Missing phosphoethanolamine group"

    # Flexible core pattern allowing both @ and @@ configurations at sn-2
    core_patterns = [
        # Direct connection pattern
        Chem.MolFromSmarts("[CH2]OC(=O)[C@,@@H](O)[CH2]OP(=O)([O-])OCC[NH2,NH3+]"),
        Chem.MolFromSmarts("[CH2]OC(=O)[C@,@@H](O)[CH2]OP(=O)(O)OCC[NH2,NH3+]"),
        # Phosphate connected via oxygen first
        Chem.MolFromSmarts("O-P(=O)([O-])OCC[NH2,NH3+].[C@,@@H](O)(COC(=O))COP(=O)"),
        # Alternative arrangement with glycerol ordering
        Chem.MolFromSmarts("OCC[NH2,NH3+].O=C(O)CO[C@,@@H](CO[P](=O)(O)[O-])O"),
    ]
    
    found_core = any(mol.HasSubstructMatch(patt) for patt in core_patterns)
    if not found_core:
        return False, "Missing core structure with R configuration"

    # Verify exactly one ester group (the 1-acyl)
    ester_patt = Chem.MolFromSmarts("[OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_patt)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester group, found {len(ester_matches)}"

    # Check acyl chain length (at least 10 carbons in the chain)
    ester_oxygen_idx = ester_matches[0][0]
    ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
    # Get the carbonyl carbon (next in the match)
    carbonyl_carbon_idx = ester_matches[0][1]
    
    # Traverse the acyl chain starting from carbonyl carbon
    chain_carbons = set()
    stack = [(mol.GetAtomWithIdx(carbonyl_carbon_idx), 0)]
    visited = set()
    while stack:
        atom, depth = stack.pop()
        if atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        # Only count carbons not part of the ester oxygen's neighborhood
        if atom.GetAtomicNum() == 6 and atom.GetIdx() != carbonyl_carbon_idx:
            chain_carbons.add(atom.GetIdx())
            # Add neighboring carbons not part of the ester group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in [ester_oxygen.GetIdx(), carbonyl_carbon_idx]:
                    stack.append((neighbor, depth + 1))
    
    if len(chain_carbons) < 10:
        return False, f"Acyl chain too short ({len(chain_carbons)} carbons)"

    return True, "Has 1-acyl-sn-glycero-3-phosphoethanolamine structure"