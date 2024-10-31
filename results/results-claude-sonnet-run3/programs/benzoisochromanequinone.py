from rdkit import Chem
from rdkit.Chem import AllChem

def is_benzoisochromanequinone(smiles: str):
    """
    Determines if a molecule is a benzoisochromanequinone.
    These are Streptomyces aromatic polyketide antibiotics with a characteristic
    benzo[g]isochromane-5,10-dione core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzoisochromanequinone, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check for required rings
        rings = mol.GetRingInfo()
        if not rings.NumRings() >= 3:
            return False, "Less than 3 rings found - needs fused ring system"

        # Look for the characteristic benzo[g]isochromane-5,10-dione pattern
        # Pattern represents fused tricyclic system with 2 carbonyls and oxygen in the right positions
        pattern = Chem.MolFromSmarts('[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6](=[O])~[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6](=[O])~[#6]~1~[#6]~O~[#6]~2')
        if not mol.HasSubstructMatch(pattern):
            return False, "Missing benzo[g]isochromane-5,10-dione core structure"

        # Additional check for aromatic character
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        if aromatic_atoms < 6:
            return False, "Insufficient aromatic character"

        # Look for carbonyl groups
        carbonyl_pattern = Chem.MolFromSmarts('C(=O)')
        if len(mol.GetSubstructMatches(carbonyl_pattern)) < 2:
            return False, "Missing required carbonyl groups"

        # Check for oxygen atom in the ring system
        oxygen_ring_pattern = Chem.MolFromSmarts('[#6]~[#8]~[#6]')
        if not mol.HasSubstructMatch(oxygen_ring_pattern):
            return False, "Missing oxygen in ring system"

        return True, "Contains benzo[g]isochromane-5,10-dione core with required features"

    except Exception as e:
        return None, f"Error in processing molecule: {str(e)}"
# Pr=None
# Recall=None