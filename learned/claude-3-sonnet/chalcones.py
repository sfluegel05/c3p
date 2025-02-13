"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:38223 chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a 1,3-diphenylpropenone (benzylideneacetophenone) or its derivatives
    formed by substitution, with the general structure Ar-CH=CH-C(=O)-Ar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for chalcone backbone pattern (Ar-CH=CH-C(=O)-Ar)
    chalcone_pattern = Chem.MolFromSmarts("[a]C=CC(=O)[a]")
    matches = mol.GetSubstructMatches(chalcone_pattern)
    if not matches:
        return False, "No chalcone backbone found"
    
    # Check for aromatic rings
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol)
                      if mol.GetRingInfo().IsAromaticRing(ring)]
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"
    
    # Check that the alpha,beta-unsaturated ketone connects the aromatic rings
    for match in matches:
        ring_atoms = set()
        for ring in aromatic_rings:
            ring_atoms.update(ring)
        if set(match).issubset(ring_atoms):
            break
    else:
        return False, "Chalcone backbone not connecting aromatic rings"
    
    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for chalcones"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 12 or o_count < 1:
        return False, "Insufficient carbon or oxygen atoms for chalcone"
    
    # Check for substituents on aromatic rings
    substituted_rings = [ring for ring in aromatic_rings
                         if any(mol.GetAtomWithIdx(idx).GetTotalNumHs() < 1
                                for idx in ring)]
    if not substituted_rings:
        return False, "No substituents found on aromatic rings"
    
    return True, "Molecule contains the chalcone backbone (Ar-CH=CH-C(=O)-Ar) with aromatic rings and substituents"