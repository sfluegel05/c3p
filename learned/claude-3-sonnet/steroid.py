"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:17240 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as 'Any of naturally occurring compounds and synthetic analogues, based on the cyclopenta[a]phenanthrene carbon skeleton, partially or completely hydrogenated; there are usually methyl groups at C-10 and C-13, and often an alkyl group at C-17. By extension, one or more bond scissions, ring expansions and/or ring contractions of the skeleton may have occurred. Natural steroids are derived biogenetically from squalene which is a triterpene.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopenta[a]phenanthrene scaffold
    scaffold = AllChem.MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_pattern = Chem.MolFromSmarts("C1CCC2C3CCCC4C5CCCC6C7CCCC%91%92%93%94%95%96%97")
    if not scaffold.HasSubstructMatch(scaffold_pattern):
        return False, "No cyclopenta[a]phenanthrene scaffold found"
    
    # Check for methyl groups at C-10 and C-13
    ring_info = mol.GetRingInfo()
    rings = [ring.AtomIds() for ring in ring_info.BondRings()]
    steroid_rings = [ring for ring in rings if len(ring) in [5, 6]]
    methyls_10_13 = False
    for ring in steroid_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_methyls = [atom for atom in ring_atoms if atom.GetTotalNumHs() == 1 and atom.GetIsAromatic() == False]
        if len(ring_methyls) == 2:
            methyls_10_13 = True
            break
    if not methyls_10_13:
        return False, "Missing methyl groups at C-10 and C-13"
    
    # Check for alkyl group at C-17 (optional)
    has_c17_alkyl = False
    for ring in steroid_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_alkyl = [atom for atom in ring_atoms if atom.GetDegree() > 3 and atom.GetIsAromatic() == False]
        if len(ring_alkyl) == 1:
            has_c17_alkyl = True
            break
    
    # Check molecular weight - steroids typically 200-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical steroid range"
    
    reason = "Contains cyclopenta[a]phenanthrene scaffold with methyl groups at C-10 and C-13"
    if has_c17_alkyl:
        reason += ", and an alkyl group at C-17"
    
    return True, reason