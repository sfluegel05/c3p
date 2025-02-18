"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Monosaccharide (parent monosaccharides are polyhydroxy aldehydes/ketones
with three or more carbons and without glycosidic connection to other sugar units)
Based on the provided definition.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A parent monosaccharide is defined as a polyhydroxy aldehyde (H[CH(OH)]nC(=O)H) or
    polyhydroxy ketone (H-[CHOH]n-C(=O)[CHOH]m-H) with three or more carbons,
    not covalently connected to another unit.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a monosaccharide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check connectivity: a monosaccharide should be a single fragment,
    # so if there are multiple disconnected parts, it is not a single unit.
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) > 1:
        return False, "Multiple fragments detected; likely not a single monosaccharide unit"
    
    # Count carbon atoms and oxygen atoms.
    atoms = mol.GetAtoms()
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}); a monosaccharide must have at least 3 carbons"
    
    # Heuristic: A general formula for many monosaccharide units is close to CnH2nOn,
    # or in the case of deoxy sugars, one oxygen less. Thus we require a minimum oxygen count
    # roughly c_count - 1 and not far above c_count. Allowing up to c_count + 2 (as in some acids).
    if not (o_count >= c_count - 1 and o_count <= c_count + 2):
        return False, f"Oxygen count ({o_count}) is not within expected range for a monosaccharide relative to {c_count} carbons"
    
    # Count hydroxyl groups using SMARTS "[OX2H]".
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(hydroxyl_smarts)
    if len(oh_matches) < 2:
        return False, f"Too few hydroxyl groups detected ({len(oh_matches)}); expect multiple -OH groups in a monosaccharide"
        
    # Look for a carbonyl group. Parent sugars carry a carbonyl (C=O) in the open-chain form.
    # However, note that cyclic sugars often do not show a free carbonyl due to hemiacetal formation.
    carbonyl_smarts = Chem.MolFromSmarts("[CX3]=[OX1]")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_smarts)
    
    # If no carbonyl is found, we assume it may be present as a latent carbonyl in a cyclic (hemiacetal) form.
    extra_note = ""
    if not has_carbonyl:
        extra_note = "No free carbonyl group detected; molecule may be in cyclic (hemiacetal) form"
    
    # Passed all heuristic tests.
    reason = (f"Matches monosaccharide criteria: {c_count} carbons, {o_count} oxygens, "
              f"{len(oh_matches)} hydroxyl groups. {extra_note}").strip()
    return True, reason
       
# Example usage:
if __name__ == "__main__":
    # Example: aldehydo-L-glucose
    test_smiles = "[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO"
    result, msg = is_monosaccharide(test_smiles)
    print(result, msg)