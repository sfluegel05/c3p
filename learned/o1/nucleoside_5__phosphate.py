"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine
    or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    
    # Purine base pattern (adenine and guanine)
    purine_smarts = 'n1c2ncnc2nc1'  # Adjusted for better accuracy
    purine_pattern = Chem.MolFromSmarts(purine_smarts)
    
    # Pyrimidine base pattern (cytosine, thymine, uracil)
    pyrimidine_smarts = 'c1c[nH]c(=O)[nH]c1=O'  # Adjusted for better accuracy
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)
    
    # Ribose sugar pattern (five-membered ring with oxygen and hydroxyl groups)
    ribose_smarts = 'O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    
    # Deoxyribose sugar pattern (ribose without hydroxyl at 2' position)
    deoxyribose_smarts = 'O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_smarts)
    
    # Phosphate group pattern attached to primary alcohol
    phosphate_smarts = 'O[P](=O)(O)O[C]'
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    
    # Validate patterns
    if None in [purine_pattern, pyrimidine_pattern, ribose_pattern, deoxyribose_pattern, phosphate_pattern]:
        return False, "Error in SMARTS pattern definitions"
    
    # Check for purine or pyrimidine base
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"
    
    # Check for ribose or deoxyribose sugar
    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not (has_ribose or has_deoxyribose):
        return False, "No ribose or deoxyribose sugar found"
    
    # Check for nucleotide linkage (base attached to sugar)
    base_match = mol.GetSubstructMatch(purine_pattern)
    if not base_match:
        base_match = mol.GetSubstructMatch(pyrimidine_pattern)
        
    sugar_match = mol.GetSubstructMatch(ribose_pattern)
    if not sugar_match:
        sugar_match = mol.GetSubstructMatch(deoxyribose_pattern)
    
    if not base_match or not sugar_match:
        return False, "Base or sugar substructure not found"
    
    # Check if base is connected to sugar via N-glycosidic bond
    connected = False
    base_atoms = set(base_match)
    sugar_atoms = set(sugar_match)
    
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (begin_idx in base_atoms and end_idx in sugar_atoms) or \
           (end_idx in base_atoms and begin_idx in sugar_atoms):
            # Check if bond is between nitrogen (base) and carbon (sugar)
            atom1 = mol.GetAtomWithIdx(begin_idx)
            atom2 = mol.GetAtomWithIdx(end_idx)
            if (atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7):
                connected = True
                break
    
    if not connected:
        return False, "Base is not connected to sugar via N-glycosidic bond"
    
    # Check for phosphate group attached to 5' carbon of sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group attached to 5' carbon of sugar"
    
    # Ensure phosphate is attached to the 5' carbon of the sugar
    phosphate_connected = False
    for match in phosphate_matches:
        phosphate_atoms = set(match)
        for atom_idx in sugar_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in phosphate_atoms:
                        phosphate_connected = True
                        break
            if phosphate_connected:
                break
        if phosphate_connected:
            break
    
    if not phosphate_connected:
        return False, "Phosphate group is not attached to 5' carbon of sugar"
    
    # At this point, criteria are satisfied
    return True, "Molecule is a nucleoside 5'-phosphate"