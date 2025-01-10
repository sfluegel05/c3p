"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for core structural elements of CoA:
    
    # 1. Adenine base - more flexible pattern
    adenine = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "Missing adenine moiety"
    
    # 2. Ribose-phosphate
    ribose_phosphate = Chem.MolFromSmarts("O[CH]1[CH][CH](O)[CH](O1)COP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose-phosphate moiety"
    
    # 3. Diphosphate
    diphosphate = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(diphosphate):
        return False, "Missing diphosphate bridge"
    
    # 4. Pantetheine with thiol
    pantetheine = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # 5. Check for 3-oxo-fatty acyl pattern
    # This is the key characteristic - a thioester connected to a beta-ketone
    oxo_pattern = Chem.MolFromSmarts("[#6]-C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo-fatty acyl pattern"

    # Get the carbon chain length
    chain_atoms = mol.GetSubstructMatch(oxo_pattern)
    if not chain_atoms:
        return False, "Cannot analyze carbon chain"
    
    # Check carbon chain characteristics
    fatty_acid_end = chain_atoms[0]  # First carbon of the pattern
    
    # Count carbons in the fatty acid chain
    def count_chain_carbons(mol, start_idx, visited=None):
        if visited is None:
            visited = set()
        visited.add(start_idx)
        count = 1
        atom = mol.GetAtomWithIdx(start_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited and neighbor.GetAtomicNum() == 6:
                count += count_chain_carbons(mol, n_idx, visited)
        return count
    
    chain_length = count_chain_carbons(mol, fatty_acid_end)
    
    if chain_length < 2:
        return False, "Fatty acid chain too short"

    # Check for presence of thioester
    thioester = Chem.MolFromSmarts("SC(=O)")
    if not mol.HasSubstructMatch(thioester):
        return False, "Missing thioester linkage"

    # Optional: Check for unsaturation (double bonds)
    double_bond = Chem.MolFromSmarts("C=C")
    has_double_bonds = mol.HasSubstructMatch(double_bond)

    # Count number of carbons and oxygens to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for 3-oxo-fatty acyl-CoA"

    return True, f"Valid 3-oxo-fatty acyl-CoA structure with {chain_length} carbons in fatty acid chain" + (" and unsaturation" if has_double_bonds else "")