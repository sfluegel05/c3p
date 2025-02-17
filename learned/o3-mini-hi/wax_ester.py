"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax ester – defined as a fatty acid ester resulting from the condensation 
of the carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.
This improved version deduplicates ester matches and ensures that the overall molecule 
contains only C and O (i.e. unfunctionalized aliphatic fragments), and that one candidate ester 
bond fragments the molecule into two pieces with long, unfunctionalized aliphatic chains.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is produced by condensation of a fatty acid (providing the carbonyl)
    with a fatty alcohol, yielding at least one ester bond that connects two long aliphatic fragments.
    The algorithm requires that the entire molecule consists only of carbon and oxygen 
    (apart from dummy atoms) and that at least one candidate ester bond (after deduplication)
    produces a fatty acid fragment (>=12 aliphatic carbons) and a fatty alcohol fragment (>= 8 aliphatic carbons).

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
    
    # Ensure the molecule contains only C and O (and dummy atoms, atomic number 0).
    for atom in mol.GetAtoms():
        atnum = atom.GetAtomicNum()
        if atnum not in (0, 6, 8):
            return False, f"Molecule contains atom with atomic number {atnum} which is not allowed for a wax ester"
            
    # Define the ester SMARTS. (This pattern may match the carbonyl and two oxygen atoms.)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    all_matches = mol.GetSubstructMatches(ester_pattern)
    if not all_matches:
        return False, "No ester group matching [CX3](=O)[OX2] found."

    # Deduplicate ester matches by collecting the actual ester bond.
    # For each match tuple (carbonyl carbon, double-bond O, linking alcohol O), retrieve the bond between the carbonyl carbon and linking oxygen.
    candidate_bond_indices = set()
    for match in all_matches:
        carbonyl_idx, dbl_ox_idx, alcohol_ox_idx = match
        bond = mol.GetBondBetweenAtoms(carbonyl_idx, alcohol_ox_idx)
        if bond is not None:
            candidate_bond_indices.add(bond.GetIdx())
            
    if not candidate_bond_indices:
        return False, "No valid ester bond found from substructure matching."
    
    # Helper: count non-aromatic aliphatic carbons in a fragment.
    def count_aliphatic_carbons(fragment):
        count = 0
        for atom in fragment.GetAtoms():
            # Only count carbon atoms, ignore dummy atoms (atomic num 0)
            if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
                count += 1
        return count

    # Helper: check if fragment contains only allowed atoms (C and O).
    def is_fragment_aliphatic(fragment):
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() not in (0, 6, 8):
                return False
        return True
    
    # Define SMARTS to identify a carbonyl group (acid fragment indicator)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)")
    
    # Minimum required aliphatic carbons in the fragments.
    min_acid_carbons = 12
    min_alcohol_carbons = 8

    # Now iterate candidate ester bond(s) to see if any candidate yields suitable fragments.
    reasons = []
    for bond_idx in candidate_bond_indices:
        # Break the molecule at this bond
        try:
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
        # Identify fragments: one should contain the carbonyl.
        for frag in frags:
            if frag.HasSubstructMatch(acid_pattern):
                frag_acid = frag
            else:
                frag_alcohol = frag
        # If not assigned by carbonyl pattern, try to assign by checking dummy atoms.
        if frag_acid is None or frag_alcohol is None:
            for frag in frags:
                has_dummy = any(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms())
                if has_dummy:
                    frag_alcohol = frag
                else:
                    frag_acid = frag
        
        if frag_acid is None or frag_alcohol is None:
            reasons.append(f"Could not clearly identify fatty acid and fatty alcohol fragments for bond {bond_idx}.")
            continue

        # Check that both fragments contain only allowed atoms.
        if not is_fragment_aliphatic(frag_acid):
            reasons.append("Fatty acid fragment contains atoms other than C and O.")
            continue
        if not is_fragment_aliphatic(frag_alcohol):
            reasons.append("Fatty alcohol fragment contains atoms other than C and O.")
            continue

        acid_carbons = count_aliphatic_carbons(frag_acid)
        alcohol_carbons = count_aliphatic_carbons(frag_alcohol)

        # Check minimum count criteria.
        if acid_carbons < min_acid_carbons:
            reasons.append(f"Fatty acid fragment too short ({acid_carbons} aliphatic C; need >= {min_acid_carbons}).")
            continue
        if alcohol_carbons < min_alcohol_carbons:
            reasons.append(f"Fatty alcohol fragment too short ({alcohol_carbons} aliphatic C; need >= {min_alcohol_carbons}).")
            continue

        # Optionally check that overall molecular weight is in the typical range for wax esters.
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 300:
            reasons.append(f"Molecular weight too low ({mol_wt:.1f} Da).")
            continue

        # If we reach here, we accept this candidate ester bond.
        return True, (f"Candidate ester bond (bond idx {bond_idx}) splits molecule into a fatty acid fragment "
                      f"with {acid_carbons} aliphatic C and a fatty alcohol fragment with {alcohol_carbons} aliphatic C.")
        
    # If we iterated through all candidate ester bonds and none fit, report the reasons.
    return False, "Wax ester not detected: " + " | ".join(reasons)

# Example usage:
# result, reason = is_wax_ester("CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC")
# print(result, reason)