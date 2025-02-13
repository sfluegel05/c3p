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
    if num_o < 8:  # Need multiple hydroxyl groups
        return False, "Insufficient oxygen atoms for hydroxyl groups"

    # Basic flavan-3-ol unit pattern (more flexible than before)
    flavan_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]=1)=[#6]-[#6](O)-[#6]=[#6]2")
    flavan_matches = len(mol.GetSubstructMatches(flavan_pattern))
    if flavan_matches < 2:
        return False, "Need at least two flavan units"

    # Check for characteristic hydroxylation patterns
    # A-ring pattern (5,7-dihydroxy)
    a_ring_pattern = Chem.MolFromSmarts("Oc1cc(O)cc2c1")
    # B-ring patterns (3',4'-dihydroxy or 3',4',5'-trihydroxy)
    b_ring_pattern1 = Chem.MolFromSmarts("c1c(O)c(O)ccc1")
    b_ring_pattern2 = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    
    if not (mol.HasSubstructMatch(a_ring_pattern) and 
            (mol.HasSubstructMatch(b_ring_pattern1) or mol.HasSubstructMatch(b_ring_pattern2))):
        return False, "Missing characteristic hydroxylation pattern"

    # Check for interflavanoid linkages
    # 4→8 linkage
    c48_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]=1)-[#6]-[#6]-[#6]-2")
    # 4→6 linkage
    c46_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]=1)-[#6]-[#6]-[#6]-2")
    
    if not (mol.HasSubstructMatch(c48_pattern) or mol.HasSubstructMatch(c46_pattern)):
        return False, "Missing characteristic interflavanoid linkage"

    # Count oxygen-containing rings (heterocycles)
    o_ring_pattern = Chem.MolFromSmarts("O1[#6][#6][#6][#6][#6]1")
    o_ring_matches = len(mol.GetSubstructMatches(o_ring_pattern))
    if o_ring_matches < 2:
        return False, "Insufficient oxygen-containing rings"

    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    if aromatic_rings < 4:  # Need at least 2 A-rings and 2 B-rings
        return False, "Insufficient aromatic rings"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Proanthocyanidins are typically large
        return False, "Molecular weight too low for proanthocyanidin"

    return True, "Contains multiple flavan units with characteristic linkages and hydroxylation patterns"