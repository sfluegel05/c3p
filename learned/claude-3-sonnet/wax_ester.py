"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: wax ester
A fatty acid ester resulting from the condensation of the carboxy group of a fatty acid 
with the alcoholic hydroxy group of a fatty alcohol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester groups found"

    # Check for elements other than C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains elements other than C, H, and O"

    # Check for aromatic rings or complex ring systems
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
        if len(ring) > 3:  # Allow only very small rings that might be computational artifacts
            return False, "Contains ring structures"

    # Look for problematic functional groups
    problematic_groups = [
        "[CX3](=O)[OX2H]",  # carboxylic acid
        "[NX3]",  # amine
        "[SX2]",  # thiol/sulfide
        "[OH1]",  # alcohol (except at chain ends)
        "O=C1OCC1",  # beta-lactone
        "O=C1OCCC1",  # gamma-lactone
        "O=C1OCCCC1",  # delta-lactone
        "c",  # aromatic carbon
        "[OX2H]C=O",  # aldehyde
        "[CX3](=O)[CX3]",  # ketone
    ]
    
    for smart in problematic_groups:
        pattern = Chem.MolFromSmarts(smart)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            # For alcohols, check if they're at chain ends
            if smart == "[OH1]":
                all_terminal = True
                for match in matches:
                    atom = mol.GetAtomWithIdx(match[0])
                    if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 1:
                        all_terminal = False
                        break
                if all_terminal:
                    continue
            return False, "Contains incompatible functional groups"

    # Check for long hydrocarbon chains on both sides of ester
    fatty_chain_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    chain_matches = mol.GetSubstructMatches(fatty_chain_pattern)
    if len(chain_matches) < 2:
        return False, "Missing long hydrocarbon chains"

    # Count carbons in main chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:  # Minimum for a very short wax ester
        return False, "Total carbon count too low for wax ester"

    # Check for excessive branching
    branching_pattern = Chem.MolFromSmarts("[CH](C)(C)C")
    if len(mol.GetSubstructMatches(branching_pattern)) > 2:
        return False, "Too highly branched for a wax ester"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    has_unsaturation = mol.HasSubstructMatch(double_bond_pattern)

    # Additional check for glycerol derivatives
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CH2X4][CHX4][CH2X4][OX2]")
    if mol.HasSubstructMatch(glycerol_pattern):
        return False, "Contains glycerol-like structure"

    return True, f"Contains ester group(s) connecting long hydrocarbon chains{' with unsaturation' if has_unsaturation else ''}"