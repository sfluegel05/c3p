"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: CHEBI:72564 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton with possible substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Core scaffold pattern for cis-fused tetrahydrochromeno[3,4-b]chromene skeleton
    # Adjusted SMARTS to match fused bicyclic system with oxygen bridges and ketone
    scaffold_pattern = Chem.MolFromSmarts("[#6]1-&@[#8]-&@[#6]-&@[#6]-&@[#6]2-&@[#6](=O)-&@[#6]-&@[#6]-&@[#6]-&@[#6]-2-&@[#6]-1")
    if not mol.HasSubstructMatch(scaffold_pattern):
        # Fallback checks for fused rings with oxygen and ketone
        # Check for ketone group
        ketone = Chem.MolFromSmarts("[CX3]=[OX1]")
        if not mol.HasSubstructMatch(ketone):
            return False, "Missing ketone group"
        
        # Check for two fused rings with oxygen atoms
        rings = mol.GetRingInfo().AtomRings()
        if len(rings) < 2:
            return False, "Insufficient rings for fused system"
        
        fused_found = False
        oxygen_in_both = False
        for i in range(len(rings)):
            for j in range(i+1, len(rings)):
                if len(set(rings[i]) & set(rings[j])) >= 2:  # Fused rings
                    fused_found = True
                    # Check oxygen in both rings
                    ring1_has_o = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in rings[i])
                    ring2_has_o = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in rings[j])
                    if ring1_has_o and ring2_has_o:
                        oxygen_in_both = True
                        break
            if oxygen_in_both:
                break
        
        if not fused_found:
            return False, "No fused rings detected"
        if not oxygen_in_both:
            return False, "Fused rings lack oxygen atoms"
        
        return True, "Fused oxygen-containing rings with ketone (fallback match)"
    
    # Verify presence of oxygen atoms in the rings
    oxygen_in_rings = any(any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms() if atom.IsInRing()) for _ in rings)
    if not oxygen_in_rings:
        return False, "Oxygen not present in ring systems"
    
    return True, "Contains core chromenochromene scaffold with ketone and oxygen bridges"