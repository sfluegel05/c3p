"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: CHEBI:66236 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids are defined as compounds consisting of a cis-fused tetrahydrochromeno[3,4-b]chromene 
    skeleton and its substituted derivatives. The term was originally restricted to natural products,
    but is now also used to describe semi-synthetic and fully synthetic compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for rotenoid core skeleton
    rotenoid_core = Chem.MolFromSmarts("[C@H]1=[C@@]2[C@@H](O[C@@H]3[C@H](c4ccc(OC)c(OC)c4O3)[C@]2(O)C(=O)c2ccccc2O1)C")
    
    # Check for rotenoid core skeleton
    if mol.HasSubstructMatch(rotenoid_core, useChirality=True):
        # Check for additional structural features
        
        # Look for cis-fused ring system
        cis_fused_rings = mol.GetRingInfo().AtomRings()
        cis_fused = any(len(ring) == 8 and Chem.AllChem.IsCisoidRing(mol, ring, 0.5) for ring in cis_fused_rings)
        if not cis_fused:
            return False, "Does not have cis-fused ring system"
        
        # Look for substituted derivatives
        substituents = [atom for atom in mol.GetAtoms() if atom.GetTotalNumHs() < atom.GetTotalValence()]
        if not substituents:
            return True, "Contains rotenoid core skeleton without additional substituents"
        else:
            return True, "Contains rotenoid core skeleton with additional substituents"
    
    # Molecule does not match rotenoid pattern
    return False, "Does not contain rotenoid core skeleton"