"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Basic molecular properties
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 30:  # Proanthocyanidins are typically large molecules
        return False, "Too small to be a proanthocyanidin"
    
    # Count key atoms
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if num_o < 8:  # Need multiple hydroxyl groups
        return False, "Insufficient oxygen atoms for hydroxyl groups"
        
    # Look for basic flavan pattern (C6-C3-C6)
    flavan_pattern = Chem.MolFromSmarts("c1cccc(CC2Oc3ccccc3CC2)c1")
    if not mol.HasSubstructMatch(flavan_pattern):
        return False, "No flavan core structure found"
    
    # Count number of aromatic rings (should have multiple)
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 4:  # Need at least 2 flavan units (4 aromatic rings)
        return False, "Insufficient aromatic rings for multiple flavan units"
    
    # Look for hydroxyl groups attached to aromatic rings
    phenol_pattern = Chem.MolFromSmarts("cO")
    phenol_matches = len(mol.GetSubstructMatches(phenol_pattern))
    
    if phenol_matches < 4:
        return False, "Insufficient phenolic hydroxyl groups"
    
    # Look for C-O-C linkages between units
    ether_pattern = Chem.MolFromSmarts("COC")
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    
    # Look for C-C linkages between units
    cc_pattern = Chem.MolFromSmarts("c1ccccc1Cc1ccccc1")
    cc_matches = len(mol.GetSubstructMatches(cc_pattern))
    
    if ether_matches < 1 and cc_matches < 1:
        return False, "No linkages found between flavan units"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Proanthocyanidins are typically large
        return False, "Molecular weight too low for proanthocyanidin"
        
    # Check for presence of multiple oxygen-containing rings
    o_ring_pattern = Chem.MolFromSmarts("O1CCCCC1")
    o_ring_matches = len(mol.GetSubstructMatches(o_ring_pattern))
    
    if o_ring_matches < 2:
        return False, "Insufficient oxygen-containing rings"
    
    return True, "Contains multiple flavan units with characteristic linkages and hydroxyl groups"