"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition: An N-glycosyl compound that has both a nucleobase 
           (normally adenine, guanine, xanthine, thymine, cytosine or uracil)
           and either a ribose or deoxyribose as functional parents.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside should include:
      - A sugar ring: typically a five-membered (furanose) ring composed of four carbons
        and one oxygen in a saturated (non-aromatic) environment.
      - A nucleobase: an aromatic heterocycle that contains at least one nitrogen.
      - An N-glycosidic bond connecting the sugar (usually at its anomeric carbon) to
        a nucleobase nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule fits the nucleoside criteria, False otherwise.
        str: A human-readable reason for the classification result.
    """
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    
    # 1. Identify the sugar ring candidate.
    # Look for a ring with exactly 5 atoms (furanose) that is non-aromatic and composed of 1 oxygen and 4 carbons.
    sugar_found = False
    sugar_indices = set()  # Will store atom indices if a sugar ring is found.
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # The sugar ring should be non-aromatic.
            if any(atom.GetIsAromatic() for atom in atoms):
                continue
            num_O = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
            num_C = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if num_O == 1 and num_C == 4:
                sugar_found = True
                sugar_indices.update(ring)
                break
    if not sugar_found:
        return False, "No ribose or deoxyribose sugar ring found (expected non-aromatic 5-membered ring with 1 oxygen and 4 carbons)"
    
    # 2. Identify the nucleobase part.
    # We expect a nucleobase to be an aromatic ring (or fused ring system) that contains at least one nitrogen.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        # Check if the ring is aromatic (all atoms in the ring are aromatic) 
        # and contains at least one nitrogen.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            num_N = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if num_N >= 1:
                nucleobase_found = True
                break
    if not nucleobase_found:
        return False, "No aromatic nucleobase ring found containing at least one nitrogen"
    
    # 3. Check for an N-glycosidic bond connecting the sugar and the nucleobase.
    # Typically, the bond is between the anomeric carbon of the sugar (which would be in the sugar ring)
    # and an aromatic nitrogen (in the nucleobase).
    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        # Check if one partner is in the sugar ring and the other is an aromatic nitrogen.
        if ((idx1 in sugar_indices and a2.GetAtomicNum() == 7 and a2.GetIsAromatic()) or 
            (idx2 in sugar_indices and a1.GetAtomicNum() == 7 and a1.GetIsAromatic())):
            glycosidic_bond_found = True
            break
    if not glycosidic_bond_found:
        return False, "No N-glycosidic bond between sugar and nucleobase detected"
    
    return True, "Molecule is a nucleoside: contains a sugar ring and an aromatic nucleobase linked via an N-glycosidic bond"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided nucleoside SMILES.
    test_smiles = "O1[C@@H](N2CN=CC=C2)[C@H](O)[C@H](O)[C@H]1CO"  # 1-(beta-D-Ribofuranosyl)-1,2-dihydropyrimidine
    result, reason = is_nucleoside(test_smiles)
    print(result, ":", reason)