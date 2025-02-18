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

    # Check for phosphoethanolamine group (O-P-O-C-C-N)
    pepatt = Chem.MolFromSmarts("[PX4]=O.OCCN")
    if not mol.HasSubstructMatch(pepatt):
        return False, "Missing phosphoethanolamine group"

    # Find glycerol backbone with correct substituents using SMARTS
    # Pattern: [sn3]-OP(=O)(O)OCCN connected to [sn2]-O and [sn1]-O-C(=O)
    # With R configuration at sn-2 (represented as @ in SMARTS)
    core_pattern = Chem.MolFromSmarts(
        "[CH2]OC(=O)"  # sn-1 acyl group
        "-[C@H](O)"    # sn-2 with R configuration
        "-[CH2]OP(=O)([O-])OCC[NH3+]"  # sn-3 phosphoethanolamine
    )
    
    if not mol.HasSubstructMatch(core_pattern):
        # Check alternative protonation states
        core_pattern2 = Chem.MolFromSmarts(
            "[CH2]OC(=O)-[C@H](O)-[CH2]OP(=O)(O)OCCN"
        )
        if not mol.HasSubstructMatch(core_pattern2):
            return False, "Missing core structure with R configuration"

    # Verify only one ester group (the 1-acyl)
    ester_patt = Chem.MolFromSmarts("[OX2]C(=O)")
    if len(mol.GetSubstructMatches(ester_patt)) != 1:
        return False, "Should have exactly one acyl group"

    # Check acyl chain length (at least 10 carbons)
    ester_carbon = mol.GetSubstructMatch(Chem.MolFromSmarts("C(=O)O"))[0]
    chain = []
    stack = [mol.GetAtomWithIdx(ester_carbon)]
    visited = set()
    
    while stack:
        atom = stack.pop()
        if atom.GetAtomicNum() == 6 and atom not in visited:
            visited.add(atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() != ester_carbon and neighbor.GetAtomicNum() == 6:
                    stack.append(neighbor)
                    chain.append(neighbor.GetIdx())

    if len(chain) < 10:
        return False, "Acyl chain too short"

    return True, "Has 1-acyl-sn-glycero-3-phosphoethanolamine structure"