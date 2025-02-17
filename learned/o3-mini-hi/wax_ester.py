"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax ester – defined as a fatty acid ester resulting from the condensation 
of the carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.
This improved version ensures the molecule contains only C and O, that exactly one candidate ester bond
splits the molecule into two acyclic (linear), unfunctionalized aliphatic fragments, 
with the fatty acid fragment having >=12 aliphatic C, and the fatty alcohol fragment having >=8 aliphatic C.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is produced by the condensation of a fatty acid (providing the carbonyl)
    with a fatty alcohol. The algorithm requires that the entire molecule consists only of carbon and oxygen,
    that at least one candidate ester bond (after deduplication) produces exactly two fragments,
    both of which are acyclic and unfunctionalized aliphatic fragments. The fragment containing a carbonyl
    must have at least 12 aliphatic carbons (fatty acid) and the other at least 8 (fatty alcohol).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a wax ester, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the parent molecule contains only allowed atoms: C (6), O (8), and optionally dummy atoms (0)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (0, 6, 8):
            return False, f"Molecule contains an atom with atomic number {atomic_num} which is not allowed for a wax ester"
    
    # Define the ester SMARTS pattern. This pattern should match the typical ester functionality.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    all_matches = mol.GetSubstructMatches(ester_pattern)
    if not all_matches:
        return False, "No ester group matching [CX3](=O)[OX2] found."
    
    # Deduplicate ester matches by retrieving the specific bond (between the carbonyl C and the alcohol O).
    candidate_bond_indices = set()
    for match in all_matches:
        # The match tuple is (carbonyl C, double-bonded O, alcohol O)
        carbonyl_idx, dbl_ox_idx, alcohol_ox_idx = match
        bond = mol.GetBondBetweenAtoms(carbonyl_idx, alcohol_ox_idx)
        if bond is not None:
            candidate_bond_indices.add(bond.GetIdx())
    
    if not candidate_bond_indices:
        return False, "No valid ester bond found from substructure matching."
    
    # Helper: count unfunctionalized aliphatic (non-aromatic) carbon atoms.
    def count_aliphatic_carbons(fragment):
        count = 0
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
                count += 1
        return count

    # Helper: Check that a fragment consists solely of C and O atoms.
    def is_fragment_aliphatic(fragment):
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() not in (0, 6, 8):
                return False
        return True

    # Helper: Check that a fragment is acyclic.
    def is_fragment_acyclic(fragment):
        ri = fragment.GetRingInfo()
        # If any rings are found, fragment is cyclic.
        return ri.NumRings() == 0

    # Define SMARTS to identify a carbonyl group so that we know which fragment is the acid part.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)")
    
    # Minimum required aliphatic carbons in the fragments.
    min_acid_carbons = 12
    min_alcohol_carbons = 8
    
    reasons = []
    # Iterate through candidate ester bonds.
    for bond_idx in candidate_bond_indices:
        try:
            # Fragment the molecule at the candidate ester bond. Dummy atoms are added.
            frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
        except Exception as e:
            reasons.append(f"Fragmentation failed at bond {bond_idx}: {str(e)}")
            continue
        
        frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        if len(frags) != 2:
            reasons.append(f"Fragmentation yielded {len(frags)} fragments (expected 2) for bond {bond_idx}.")
            continue
        
        frag_acid = None
        frag_alcohol = None
        # Identify fragments: the one that contains the carbonyl is the acid fragment.
        for frag in frags:
            if frag.HasSubstructMatch(acid_pattern):
                frag_acid = frag
            else:
                frag_alcohol = frag
        # If not clearly identified, use dummy atom presence as a secondary means.
        if frag_acid is None or frag_alcohol is None:
            for frag in frags:
                if any(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms()):
                    frag_alcohol = frag
                else:
                    frag_acid = frag
        
        if frag_acid is None or frag_alcohol is None:
            reasons.append(f"Could not clearly identify fatty acid and fatty alcohol fragments for bond {bond_idx}.")
            continue

        # Both fragments must contain only C and O.
        if not is_fragment_aliphatic(frag_acid):
            reasons.append("Fatty acid fragment contains atoms other than C and O.")
            continue
        if not is_fragment_aliphatic(frag_alcohol):
            reasons.append("Fatty alcohol fragment contains atoms other than C and O.")
            continue
        
        # Both fragments should be acyclic (no rings).
        if not is_fragment_acyclic(frag_acid):
            reasons.append("Fatty acid fragment is cyclic; expected a linear aliphatic chain.")
            continue
        if not is_fragment_acyclic(frag_alcohol):
            reasons.append("Fatty alcohol fragment is cyclic; expected a linear aliphatic chain.")
            continue

        acid_carbons = count_aliphatic_carbons(frag_acid)
        alcohol_carbons = count_aliphatic_carbons(frag_alcohol)
        
        if acid_carbons < min_acid_carbons:
            reasons.append(f"Fatty acid fragment too short ({acid_carbons} aliphatic C; need >= {min_acid_carbons}).")
            continue
        if alcohol_carbons < min_alcohol_carbons:
            reasons.append(f"Fatty alcohol fragment too short ({alcohol_carbons} aliphatic C; need >= {min_alcohol_carbons}).")
            continue
        
        # Optional: Check that the overall molecular weight meets a minimum size criterion.
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 300:
            reasons.append(f"Molecular weight too low ({mol_wt:.1f} Da).")
            continue
        
        # If we reach here, we have found a candidate ester bond that splits the molecule appropriately.
        return True, (f"Candidate ester bond (bond idx {bond_idx}) splits molecule into a fatty acid fragment "
                      f"with {acid_carbons} aliphatic C and a fatty alcohol fragment with {alcohol_carbons} aliphatic C.")
    
    return False, "Wax ester not detected: " + " | ".join(reasons)

# Example usage:
# result, reason = is_wax_ester("CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC")
# print(result, reason)