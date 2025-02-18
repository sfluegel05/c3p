"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a hydroxyflavan with a flavan-3-ol skeleton and substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for flavan-3-ol core: oxygen in a six-membered ring with adjacent hydroxyl
    core_pattern = Chem.MolFromSmarts("[O;R]C[C;R](O)")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No flavan-3-ol core (O-C-C-OH in ring) found"
    
    # Check for two aromatic rings (A and B rings)
    aromatic_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) >= 6 and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            aromatic_rings.append(ring)
    if len(aromatic_rings) < 2:
        return False, f"Found {len(aromatic_rings)} aromatic rings, need at least 2"
    
    # Verify hydroxyl count (at least 3 for hydroxyflavan class)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 3:
        return False, f"Only {hydroxyl_count} hydroxyl groups found"
    
    # Check molecular weight range (typical catechins 290-450 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical catechin range (200-1000 Da)"
    
    return True, "Contains flavan-3-ol core with two aromatic rings and multiple hydroxyl groups"