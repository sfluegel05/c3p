"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid

Definition:
  Any dicarboxylic acid carrying a hydroxy group on the carbon atom
  at position alpha to a carboxyl group.

Algorithm improvements over the previous approach:
  1. Adds explicit hydrogens (to help distinguish –OH from carbonyl oxygens).
  2. Requires exactly 2 carboxyl groups (using the SMARTS "[CX3](=O)[OX2H]").
  3. For each acid, inspects directly bonded carbons (alpha candidates) that are
     not in rings and not aromatic.
  4. Requires the candidate alpha carbon to have at least one genuine hydroxyl 
     substituent (an oxygen atom with at least one hydrogen via a single bond)
     and to not also have a double-bonded (carbonyl) oxygen aside from the carboxyl attachment.
  5. Rejects molecules with nitrogen atoms.
"""

from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    
    A 2-hydroxydicarboxylic acid is defined as any dicarboxylic acid 
    carrying a hydroxy group on the carbon atom at position alpha to a carboxyl group.
    
    The algorithm:
      1. Parses and adds explicit hydrogens.
      2. Rejects molecules with nitrogen atoms.
      3. Searches for exactly two carboxyl groups (using SMARTS).
      4. For each carboxyl group, finds a candidate alpha-carbon (a carbon bonded to the acid carbon,
         not aromatic and not in a ring).
      5. Checks that this candidate bears a hydroxyl substituent (an O with at least one H attached)
         and that the candidate does not carry any additional carbonyl (C=O) oxygen.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies, False otherwise.
        str : Explanation for the decision.
    """
    # Parse SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Reject molecules containing nitrogen.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains unexpected nitrogen atoms"
    
    # SMARTS pattern for carboxyl group (–COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    if carboxyl_pattern is None:
        return False, "Error generating carboxyl SMARTS pattern"
    
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl groups found"
    # We take the acid carbon as the first atom in each match.
    acid_carbon_indices = set(match[0] for match in carboxyl_matches)
    if len(acid_carbon_indices) != 2:
        return False, f"Molecule has {len(acid_carbon_indices)} carboxyl groups (exactly 2 required)"
    
    # Helper function: determine if an oxygen atom is a genuine hydroxyl (has at least one hydrogen attached)
    def is_hydroxyl(oxygen_atom):
        if oxygen_atom.GetAtomicNum() != 8:
            return False
        # Check for at least one hydrogen neighbor:
        return any(neigh.GetAtomicNum() == 1 for neigh in oxygen_atom.GetNeighbors())
    
    # For each acid carbon, check its neighbors to find an alpha candidate.
    for acid_idx in acid_carbon_indices:
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # Inspect neighbors of the acid carbon.
        for neighbor in acid_atom.GetNeighbors():
            # Consider only carbon atoms that are not themselves acid carbons.
            if neighbor.GetAtomicNum() != 6 or neighbor.GetIdx() in acid_carbon_indices:
                continue
            # Require candidate alpha carbon is not in a ring and not aromatic.
            if neighbor.IsInRing() or neighbor.GetIsAromatic():
                continue
            
            # Check that the candidate alpha carbon has at least one hydroxyl substituent.
            hydroxyl_found = False
            for subnbr in neighbor.GetNeighbors():
                # Skip the acid carbon connection.
                if subnbr.GetIdx() == acid_idx:
                    continue
                # Check if this neighbor is oxygen and is a hydroxyl.
                if subnbr.GetAtomicNum() == 8 and is_hydroxyl(subnbr):
                    hydroxyl_found = True
                    break
            if not hydroxyl_found:
                continue
            
            # Also verify that the candidate alpha carbon does not have any extra carbonyl (C=O) oxygen.
            # (We look at bonds from the candidate and if any oxygen is double-bonded, we discount it.)
            extra_carbonyl = False
            for bond in neighbor.GetBonds():
                # Skip bond back to acid carbon (our intended connectivity).
                other = bond.GetOtherAtom(neighbor)
                if other.GetIdx() == acid_idx:
                    continue
                if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    extra_carbonyl = True
                    break
            if extra_carbonyl:
                continue
            
            # If we reach here, we found a candidate alpha-carbon that is attached to a hydroxyl and
            # does not bear an extra carbonyl group.
            return True, "Molecule is a 2-hydroxydicarboxylic acid with an alpha hydroxy substituent"
    
    return False, "No suitable alpha-carbon with a free hydroxyl substituent adjacent to a carboxyl group was found"


# Uncomment below to run tests
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
        ("OC(C(C1CC1=C)C(O)=O)C(O)=O", "2-hydroxy-3-(2-methylidenecyclopropyl)butanedioic acid"),
        ("CC[C@@](O)(CC(O)=O)C(O)=O", "(R)-2-ethylmalic acid"),
        ("OC(CCCCCC(O)=O)C(O)=O", "2-hydroxyoctanedioic acid"),
        ("C(C(CP(=O)(O)[H])(C(O)=O)O)C(O)=O", "2-phosphinomethylmalic acid"),
        ("CC(C(O)=O)C(C)(O)C(O)=O", "2,3-dimethylmalic acid"),
        ("CC(C)C(C(O)C(O)=O)C(O)=O", "3-isopropylmalic acid"),
        ("C[C@](O)(CC(O)=O)C(O)=O", "L-citramalic acid"),
        ("O[C@@H](CC(=O)C(O)=O)[C@@H](O)C(O)=O", "5-dehydro-4-deoxy-D-glucaric acid"),
        ("C(C(CP(O)=O)C(O)=O)(C(O)=O)O", "phosphinomethylisomalic acid"),
        ("OC(=O)\\C=C/C=C(/O)C(O)=O", "(2E,4Z)-2-hydroxymuconic acid"),
        ("OC(=O)\\C=C\\C=C(/O)C(O)=O", "(2Z,4E)-2-hydroxymuconic acid"),
        ("CC(O)(CC(O)=O)C(O)=O", "citramalic acid"),
        ("C(C(C(O)=O)O)(CCSC)C(=O)O", "3-(2-methylthioethyl)malic acid"),
        ("OC(CCCC(O)=O)C(O)=O", "2-hydroxyadipic acid"),
        # False positives:
        ("OC(CC1=CC=C(O)C=C1)(C(O)C(O)=O)C(O)=O", "piscidic acid"),
        ("CSCCCCCC(C(O)C(O)=O)C(O)=O", "3-(5'-Methylthio)pentylmalic acid"),
        ("O[C@H](CC(=O)C(O)=O)C(O)=O", "D-4-hydroxy-2-oxoglutaric acid"),
    ]
    
    for smi, name in test_examples:
        res, reason = is_2_hydroxydicarboxylic_acid(smi)
        print(f"{name}: {res}, Reason: {reason}")