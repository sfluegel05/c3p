"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid
Definition: Any dicarboxylic acid carrying a hydroxy group on the carbon atom 
at position alpha to a carboxy group.
The algorithm requires that the molecule:
  1. Has exactly 2 carboxyl (–COOH) groups, identified with the SMARTS pattern "[CX3](=O)[OX2H]".
  2. Contains no nitrogen atoms (to avoid peptides and decorated molecules).
  3. Has at least one carboxyl group whose carbon (the –COOH carbon) is directly 
     connected to a saturated (sp3) carbon that carries an –OH (with at least one attached hydrogen).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise.
        str : Reason for classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help with substructure matching.
    mol = Chem.AddHs(mol)
    
    # Reject any molecule that contains nitrogen, since no true positives have N.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains nitrogen atoms, which is not expected for a 2-hydroxydicarboxylic acid"
    
    # Define a SMARTS pattern for carboxyl (–COOH) groups.
    # The pattern "[CX3](=O)[OX2H]" matches a carbonyl carbon with a hydroxyl oxygen carrying hydrogen.
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    if carboxyl_pattern is None:
        return False, "Error creating carboxyl SMARTS pattern"
    
    # Find substructure matches for carboxyl groups.
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl groups found"
    
    # Deduplicate matches by the carboxyl carbon index (first atom in the SMARTS pattern).
    carboxyl_carbon_idxs = set(match[0] for match in carboxyl_matches)
    if len(carboxyl_carbon_idxs) != 2:
        return False, f"Molecule has {len(carboxyl_carbon_idxs)} carboxyl groups (exactly 2 required)"
    
    # Define a helper function to check if an atom is a proper hydroxyl oxygen.
    def is_hydroxyl(oxygen_atom):
        # Should be oxygen and have at least one hydrogen as neighbor.
        if oxygen_atom.GetAtomicNum() != 8:
            return False
        return any(neigh.GetAtomicNum() == 1 for neigh in oxygen_atom.GetNeighbors())
    
    # Look for an alpha (neighboring) carbon to any carboxyl carbon that is sp3 and bears an –OH.
    alpha_found = False
    for acid_carbon_idx in carboxyl_carbon_idxs:
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        # Iterate over neighbors; we only consider carbon atoms.
        for neighbor in acid_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            # Check that the candidate alpha carbon is sp3.
            if neighbor.GetHybridization() != rdchem.HybridizationType.SP3:
                continue
            # Check if the candidate alpha carbon carries a hydroxyl group (other than the acid's oxygens)
            for subnbr in neighbor.GetNeighbors():
                # Skip if it is the acid carbon itself.
                if subnbr.GetIdx() == acid_carbon_idx:
                    continue
                if is_hydroxyl(subnbr):
                    alpha_found = True
                    break
            if alpha_found:
                break
        if alpha_found:
            break

    if not alpha_found:
        return False, "No alpha carbon adjacent to a carboxyl group carries an sp3 hydroxyl substituent"
    
    return True, "Molecule is a 2-hydroxydicarboxylic acid with an alpha hydroxy substituent"

# Uncomment the following test block to run some examples.
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("C[C@H](C(O)=O)[C@@](C)(O)C(O)=O", "(2R,3S)-2,3-dimethylmalic acid"),
        ("O[C@@H](CCC(O)=O)C(O)=O", "(S)-2-hydroxyglutaric acid"),
        ("O[C@H](CC(O)=O)C(O)=O", "(R)-malic acid"),
        ("CCC(C(O)C(O)=O)C(O)=O", "3-ethylmalic acid"),
        ("CC(C)[C@@H]([C@@H](O)C(O)=O)C(O)=O", "(2R,3S)-3-isopropylmalic acid"),
        ("CC(C)([C@@H](O)C(O)=O)C(O)=O", "(R)-3,3-dimethylmalic acid"),
        ("C[C@@](O)(CC(O)=O)C(O)=O", "D-citramalic acid"),
        ("OC(CCC(O)=O)C(O)=O", "2-hydroxyglutaric acid"),
        # False positives (should return False):
        ("OC(CC(N)C(O)=O)C(O)=O", "(4R)-4-Hydroxy-L-glutamic acid"),
        ("O[C@H]1C=C(C=C[C@]1(O)C(O)=O)C(O)=O", "(3S,4R)-3,4-dihydroxycyclohexa-1,5-diene-1,4-dicarboxylic acid"),
    ]
    for smi, name in test_examples:
        result, reason = is_2_hydroxydicarboxylic_acid(smi)
        print(f"{name}: {result}, Reason: {reason}")