"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: alkaloids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms present"
        
    # Check for heterocyclic rings containing nitrogen
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings present"
        
    # Find rings containing nitrogen
    has_n_ring = False
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            has_n_ring = True
            break
            
    if not has_n_ring:
        return False, "No nitrogen-containing rings found"
    
    # Check for peptide/protein characteristics
    peptide_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(peptide_pattern):
        peptide_matches = len(mol.GetSubstructMatches(peptide_pattern))
        if peptide_matches > 1:
            return False, "Appears to be a peptide"
    
    # Check for nucleotide characteristics
    nucleotide_pattern = Chem.MolFromSmarts("[$(O[P](=O)(O)O)]")
    if mol.HasSubstructMatch(nucleotide_pattern):
        return False, "Contains phosphate group characteristic of nucleotides"
        
    # Check for basic nitrogen (not amide)
    basic_n_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N~P)]")
    if not mol.HasSubstructMatch(basic_n_pattern):
        return False, "No basic nitrogen atoms found"
        
    # Natural product-like characteristics
    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical alkaloid range"
        
    # Check for reasonable number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1 or num_rings > 10:
        return False, f"Number of rings ({num_rings}) outside typical alkaloid range"
        
    # Check for reasonable number of stereocenters (common in natural products)
    num_stereo = len(Chem.FindMolChiralCenters(mol))
    
    # Calculate fraction of sp3 carbons (natural products tend to be more sp3-rich)
    num_sp3 = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    sp3_fraction = num_sp3/total_c if total_c > 0 else 0
    
    # Natural products tend to have more complex ring systems and stereocenters
    if num_stereo == 0 and sp3_fraction < 0.2:
        return False, "Structure lacks complexity typical of alkaloids"
        
    return True, "Contains nitrogen heterocycle(s) with basic nitrogen and natural product-like features"