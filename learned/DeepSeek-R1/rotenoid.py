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
    
    # Check for essential ketone group
    ketone = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(ketone):
        return False, "Missing ketone group"
    
    # Get all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    
    # Check for core scaffold: cis-fused tetrahydrochromeno[3,4-b]chromene
    # SMARTS pattern attempts to capture fused bicyclic system with oxygen bridges
    scaffold_pattern = Chem.MolFromSmarts("[#6]1-&@[#8]-&@[#6]-&@[#6]-&@[#6]2-&@[#6](=O)-&@[#6]-&@[#6]-&@[#6]-&@[#6]-2-&@[#6]-1")
    if mol.HasSubstructMatch(scaffold_pattern):
        # Verify oxygen presence in the matched scaffold
        matched_atoms = {atom for match in mol.GetSubstructMatches(scaffold_pattern) for atom in match}
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in matched_atoms):
            return True, "Core scaffold matched with oxygen and ketone"
    
    # Fallback check for fused oxygen-containing rings
    fused_oxygen_found = False
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # Check for fused rings (sharing at least 2 atoms)
            if len(set(rings[i]) & set(rings[j])) >= 2:
                # Check both rings contain oxygen
                ring1_o = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in rings[i])
                ring2_o = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in rings[j])
                if ring1_o and ring2_o:
                    fused_oxygen_found = True
                    break
        if fused_oxygen_found:
            break
    
    if not fused_oxygen_found:
        return False, "No fused oxygen-containing rings"
    
    return True, "Fused oxygen-containing rings with ketone group"