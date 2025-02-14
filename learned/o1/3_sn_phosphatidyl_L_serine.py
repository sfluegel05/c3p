"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:64381 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine has a glycerol backbone with acyl groups at positions 1 and 2,
    and a phospho-L-serine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define phospho-L-serine pattern
    phosphoserine_smarts = '[O]-P(=O)([O-])-[O]-C[C@@H](N)C(=O)O'
    phosphoserine_pattern = Chem.MolFromSmarts(phosphoserine_smarts)
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Phospho-L-serine group not found"

    # Define glycerol backbone pattern (three carbons each attached to oxygen)
    glycerol_smarts = '[CH2]-[CH]-[CH2]'
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "Glycerol backbone not found"

    # Check for ester linkages at positions sn-1 and sn-2
    ester_pattern = Chem.MolFromSmarts('[C;H2](OC(=O)[C])[C;H](OC(=O)[C])[CH2]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "Ester linkages at sn-1 and sn-2 not found"

    # Verify that glycerol is connected to phospho-L-serine at sn-3 position
    # Define the linkage between glycerol backbone and phospho-L-serine
    glycerol_phosphate_smarts = '[C;H2][C;H](O[P](=O)(O)[O][CH2][C@@H](N)C(=O)O)[CH2]'
    glycerol_phosphate_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol backbone not connected to phospho-L-serine at sn-3"

    # Check acyl chain lengths to confirm they are fatty acids
    acyl_chain_lengths = []
    ester_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)O[CH]'))
    if len(ester_bonds) < 2:
        return False, "Less than two acyl chains found"
    for bond in ester_bonds:
        acyl_carbon_idx = bond[0]  # Carbonyl carbon
        # Trace the chain from the carbonyl carbon
        chain_length = 0
        visited = set()
        stack = [acyl_carbon_idx]
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                for neighbor in atom.GetNeighbors():
                    nbr_idx = neighbor.GetIdx()
                    nbr_atomic_num = neighbor.GetAtomicNum()
                    if nbr_idx not in visited and nbr_atomic_num in (6, 1):  # Continue through carbons and hydrogens
                        stack.append(nbr_idx)
        acyl_chain_lengths.append(chain_length)

    if min(acyl_chain_lengths) < 8:
        return False, "Acyl chains are too short to be fatty acids"

    return True, "Molecule is a 3-sn-phosphatidyl-L-serine with correct glycerol backbone and substituents"