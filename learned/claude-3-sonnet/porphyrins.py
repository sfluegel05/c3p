"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: porphyrins
Natural pigments containing a fundamental skeleton of four pyrrole nuclei united through 
the alpha-positions by four methine groups to form a macrocyclic structure.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_porphyrin, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for the basic porphyrin core pattern:
    # Four pyrrole rings connected by methine bridges
    porphyrin_core = Chem.MolFromSmarts("[nH,n-,n+0]1ccc2c1[c,C]=c1[c,C]=c3[nH,n-,n+0]c(cc4[nH,n-,n+0]c(cc5[nH,n-,n+0]c(c1)cc5)cc4)cc3c2")
    metal_porphyrin_core = Chem.MolFromSmarts("[nX2,nX2-,nX2+0]1ccc2c1[c,C]=c1[c,C]=c3[nX2,nX2-,nX2+0]c(cc4[nX2,nX2-,nX2+0]c(cc5[nX2,nX2-,nX2+0]c(c1)cc5)cc4)cc3c2")
    
    has_core = mol.HasSubstructMatch(porphyrin_core) or mol.HasSubstructMatch(metal_porphyrin_core)
    if not has_core:
        return False, "No porphyrin core structure found"

    # Count nitrogens - should have at least 4 for the core structure
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient nitrogen atoms for porphyrin structure"
        
    # Look for aromatic system
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 20:  # Core should have at least 20 aromatic atoms
        return False, "Insufficient aromatic system for porphyrin"
        
    # Check for presence of metal coordination (optional)
    metals = ["Fe", "Mg", "Zn", "Co", "Ni", "Cu", "Pd", "Pt"]
    has_metal = any(atom.GetSymbol() in metals for atom in mol.GetAtoms())
    
    # Verify macrocyclic structure by checking ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 5:  # Should have at least 5 rings (4 pyrroles + central)
        return False, "Insufficient ring count for porphyrin structure"
        
    # Additional check for conjugated system
    conjugated_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Lacks required conjugated system"

    if has_metal:
        return True, "Metalloporphyrin with characteristic tetrapyrrole macrocycle"
    else:
        return True, "Free-base porphyrin with characteristic tetrapyrrole macrocycle"