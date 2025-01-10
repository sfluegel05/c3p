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

    # Normalize charges - convert to neutral form where possible
    mol = Chem.AddHs(mol)
    
    # Check for CoA core structure
    coA_pattern = Chem.MolFromSmarts("[$(N1C=NC2=C1N=CN=C2N)]" + # Adenine
                                    "[$(C[C@H]1O[C@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(O)(O)=O)]" + # Ribose phosphate
                                    "[$(CC(C)(C)COP(O)(=O)OP(O)(=O)O)]" + # Diphosphate
                                    "[$(SCCNC(=O)CCNC(=O))]") # Pantetheine
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Missing or incorrect CoA structure"

    # Check for specific 3-oxo-fatty acyl pattern
    # This enforces the exact position of the oxo group and linear nature
    oxo_pattern = Chem.MolFromSmarts("[#6]-C(=O)CC(=O)SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing or incorrect 3-oxo-fatty acyl pattern"

    # Exclude molecules with cyclic structures in the fatty acid portion
    # First, find the thioester carbon
    thioester = Chem.MolFromSmarts("C(=O)SCCNC")
    matches = mol.GetSubstructMatches(thioester)
    if not matches:
        return False, "Cannot locate thioester group"
    
    # Get the carbon atom index of the thioester
    thioester_idx = matches[0][0]
    
    # Check for rings that include atoms before the thioester
    ri = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ri.AtomRings():
        ring_atoms.update(ring)
    
    # If any ring atoms are connected to the fatty acid portion, reject
    for atom_idx in ring_atoms:
        if atom_idx < thioester_idx:
            return False, "Contains cyclic structures in fatty acid portion"

    # Check for branching in the fatty acid portion
    branched_pattern = Chem.MolFromSmarts("[CH2][CH]([#6])[#6]-C(=O)CC(=O)S")
    if mol.HasSubstructMatch(branched_pattern):
        return False, "Contains branched structures in fatty acid portion"

    # Verify linear chain length
    chain_pattern = Chem.MolFromSmarts("CCCC-C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Fatty acid chain too short"

    # Additional check for aromatic rings in the fatty acid portion
    aromatic_pattern = Chem.MolFromSmarts("a-C(=O)CC(=O)S")
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings in fatty acid portion"

    # Check molecular weight is in reasonable range for 3-oxo-fatty acyl-CoA
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if not (850 < mol_wt < 1200):
        return False, f"Molecular weight {mol_wt:.1f} outside expected range for 3-oxo-fatty acyl-CoA"

    return True, "Valid 3-oxo-fatty acyl-CoA structure with correct CoA moiety, 3-oxo group, and linear fatty acid chain"