from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vinca_alkaloid(smiles: str):
    """
    Determines if a molecule is a vinca alkaloid based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a vinca alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for presence of indole core
    indole_pattern = Chem.MolFromSmarts('[#6]1:[#6]:[#6]:[#6]:[#6]2:[#6]:1:[nH]:[#6]:[#6]:2')
    
    # Check for presence of indoline core (reduced indole)
    indoline_pattern = Chem.MolFromSmarts('[#6]1:[#6]:[#6]:[#6]:[#6]2:[#6]:1:[N]:[#6]:[#6]:2')
    
    # Need at least one of indole or indoline
    has_indole = mol.HasSubstructMatch(indole_pattern)
    has_indoline = mol.HasSubstructMatch(indoline_pattern)
    
    if not (has_indole or has_indoline):
        return False, "Missing required indole/indoline core structure"

    # Check for presence of bridged ring system
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 5:
        return False, "Insufficient ring systems for vinca alkaloid"
        
    # Check for presence of ester group (common in vinca alkaloids)
    ester_pattern = Chem.MolFromSmarts('C(=O)OC')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing characteristic ester group"
        
    # Check molecular weight range typical for vinca alkaloids
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if not (400 < mol_weight < 1200):
        return False, f"Molecular weight {mol_weight:.1f} outside typical range for vinca alkaloids"

    # Check for nitrogen content (alkaloids should have multiple N atoms)
    num_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_nitrogens < 2:
        return False, "Insufficient nitrogen atoms for vinca alkaloid structure"
        
    # Check for characteristic dimeric structure by looking for two nitrogen-containing ring systems
    # separated by appropriate distance
    n_rings = mol.GetRingInfo().AtomRings()
    n_containing_rings = []
    
    for ring in n_rings:
        for atom_idx in ring:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7:
                n_containing_rings.append(ring)
                break
                
    if len(n_containing_rings) < 2:
        return False, "Missing characteristic dimeric ring structure"
        
    return True, "Structure contains characteristic vinca alkaloid features"
# Pr=1.0
# Recall=0.6956521739130435