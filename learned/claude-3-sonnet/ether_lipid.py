"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:38116 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a lipid with one or more ether linkages between the glycerol backbone and alkyl chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern ([CH2X4][CHX4][CH2X4] with 2 or 3 ether/ester oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    # Check if at least one glycerol match has 2 or 3 ether/ester oxygens
    for match in glycerol_matches:
        glycerol_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        ether_ester_oxygens = sum(1 for atom in glycerol_atoms if atom.GetTotalDegree() == 2 and atom.GetAtomicNum() == 8)
        if ether_ester_oxygens in [2, 3]:
            break  # Found a valid glycerol backbone
    else:
        return False, "No valid glycerol backbone found"

    # Look for ether linkages (-O-C)
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages found"

    # Look for long alkyl chains (>= 6 carbons)
    alkyl_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    if not alkyl_matches:
        return False, "No long alkyl chains found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short, need at least 8 rotatable bonds"

    # Check molecular weight - ether lipids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for ether lipid"

    return True, "Contains glycerol backbone with one or more ether-linked alkyl chains"