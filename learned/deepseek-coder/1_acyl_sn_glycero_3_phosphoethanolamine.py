"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: CHEBI:17504 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoethanolamine has a glycerol backbone with (R)-configuration,
    a single acyl chain at the 1-position, and a phosphoethanolamine group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core pattern with strict stereochemistry
    core_pattern = Chem.MolFromSmarts("[C@@H](O)(COC(=O)[*])COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No glycerol backbone with (R)-configuration found"

    # Count ester bonds - should be exactly 1 for monoacyl
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester bonds, need exactly 1"

    # Verify single acyl chain by checking for only one carbonyl group connected to oxygen
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl chains, need exactly 1"

    # Check acyl chain length - should be at least 10 carbons
    acyl_chain = mol.GetSubstructMatch(acyl_pattern)
    if acyl_chain:
        acyl_atom = mol.GetAtomWithIdx(acyl_chain[0])
        chain_length = 0
        for neighbor in acyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                chain_length = max(chain_length, len(Chem.rdmolops.GetShortestPath(mol, neighbor.GetIdx(), -1)))
        if chain_length < 10:
            return False, f"Acyl chain too short: {chain_length} carbons, need at least 10"

    # Check molecular weight - typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for 1-acyl-sn-glycero-3-phosphoethanolamine"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for 1-acyl-sn-glycero-3-phosphoethanolamine"
    if o_count != 6:
        return False, "Must have exactly 6 oxygens (1 ester, 1 phosphate, 2 hydroxyls, 1 ether, 1 amine)"

    return True, "Contains glycerol backbone with (R)-configuration, single acyl chain at 1-position, and phosphoethanolamine group at 3-position"