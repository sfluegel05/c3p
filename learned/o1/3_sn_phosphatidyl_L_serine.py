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

    # Define phospho-L-serine pattern (adjusted without charges and specific stereochemistry)
    phosphoserine_smarts = '[O]-P(=O)(O)-[O]-C-C(N)C(=O)O'
    phosphoserine_pattern = Chem.MolFromSmarts(phosphoserine_smarts)
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Phospho-L-serine group not found"

    # Define glycerol backbone pattern (three carbons connected in a chain)
    glycerol_smarts = '[CH2]-[CH]-[CH2]'
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "Glycerol backbone not found"

    # Check for ester linkages at positions sn-1 and sn-2
    # Ester linkage pattern: [O]-C(=O)-[C]
    ester_smarts = '[O]-C(=O)-[C]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Less than two ester linkages found"

    # Verify that glycerol is connected to phospho-L-serine at sn-3 position
    # Define the linkage between glycerol backbone and phospho-L-serine
    glycerol_phosphate_smarts = '[CH2][CH](O[P](=O)(O)O[CH2][CH](N)C(=O)O)[CH2]'
    glycerol_phosphate_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol backbone not connected to phospho-L-serine at sn-3"

    # Check acyl chain lengths to confirm they are fatty acids
    # Identify acyl chains attached via ester linkages
    acyl_chain_lengths = []
    for match in ester_matches:
        ester_oxygen_idx = match[0]  # Index of the oxygen atom in the ester linkage
        carbonyl_carbon_idx = mol.GetAtomWithIdx(match[1]).GetIdx()  # Carbonyl carbon
        # Walk the chain from the carbonyl carbon
        chain_length = 0
        visited = set()
        stack = [carbonyl_carbon_idx]
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