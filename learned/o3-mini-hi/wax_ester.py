"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax ester – defined as a fatty acid ester resulting from the condensation 
of the carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

This improved version requires:
  • The molecule contains only allowed atoms (C and O).
  • Exactly one ester group is present.
  • Fragmentation at the ester bond yields exactly two acyclic fragments.
  • The fragment containing the carbonyl (acid fragment) has a linear, unfunctionalized
    aliphatic chain with >= 12 carbons.
  • The other fragment (alcohol fragment) is a linear, unfunctionalized aliphatic chain with >= 8 carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    
    A wax ester is produced by the condensation of a fatty acid and a fatty alcohol.
    This function requires that the molecule (i) contains only carbon and oxygen,
    (ii) exactly one ester group ([CX3](=O)[OX2]) is present, and (iii) fragmentation
    at the ester bond produces two acyclic fragments. One of these fragments must contain 
    the carbonyl and represent a fatty acid chain (>=12 aliphatic carbons) while the 
    other is a fatty alcohol chain (>=8 aliphatic carbons). Moreover, the fragments must be 
    “unfunctionalized” (i.e. linear with minimal branching).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a wax ester, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check that the molecule contains only allowed atoms: C (6) and O (8)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (0, 6, 8):  # dummy atoms (0) allowed from fragmentation
            return False, f"Disallowed atom with atomic number {atom.GetAtomicNum()} detected."
            
    # Define the ester SMARTS pattern.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester group ([CX3](=O)[OX2]) found."
    
    # Require that there is exactly one ester group.
    if len(matches) != 1:
        return False, f"Expected exactly one ester group but found {len(matches)} matches."
    
    # Deduplicate candidate ester bonds by retrieving the bond connecting the carbonyl C and the alcohol O.
    candidate_bond_indices = set()
    # For our pattern the match tuple is (carbonyl C, carbonyl O (double-bonded), alcohol O)
    for match in matches:
        carbonyl_idx, _, alcohol_ox_idx = match
        bond = mol.GetBondBetweenAtoms(carbonyl_idx, alcohol_ox_idx)
        if bond is not None:
            candidate_bond_indices.add(bond.GetIdx())
    if len(candidate_bond_indices) != 1:
        return False, f"Expected exactly one valid ester bond from substructure matching, found {len(candidate_bond_indices)}."
    
    bond_idx = candidate_bond_indices.pop()
    
    # Helper functions
    
    def count_aliphatic_carbons(fragment):
        """Count non-aromatic carbon atoms (ignoring dummy atoms)."""
        return sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic())
    
    def is_fragment_aliphatic(fragment):
        """Ensure fragment atoms (ignoring dummy atoms) are only C or O."""
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() == 0:
                continue
            if atom.GetAtomicNum() not in (6, 8):
                return False
        return True
    
    def is_fragment_acyclic(fragment):
        """Check that a fragment is acyclic."""
        return fragment.GetRingInfo().NumRings() == 0
    
    def remove_dummies(fragment):
        """Return list of atoms (and bonds) ignoring dummy (atomic num 0) atoms.
           We do not reconstruct a new molecule, but use this for connectivity checks.
        """
        return [atom for atom in fragment.GetAtoms() if atom.GetAtomicNum() != 0]
    
    def is_linear_aliphatic(fragment, role):
        """
        Check that the fragment is a linear (unbranched) aliphatic chain.
        For each carbon (ignoring dummy atoms), count connections to other carbons.
        A non-terminal carbon should have exactly 2 carbon neighbors.
        For role 'acid', we allow one carbonyl oxygen (which should be connected by a double bond).
        For role 'alcohol', any O present must be terminal (i.e. attached to exactly one carbon).
        """
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() == 0:
                continue
            if atom.GetAtomicNum() == 6:  # carbon
                # Count neighbors that are carbons (ignore dummy atoms)
                c_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
                # Terminal carbon will have 1 C-neighbor, internal should have 2.
                if len(c_neighbors) > 2:
                    return False
            elif atom.GetAtomicNum() == 8:
                # For acid fragment: allow if oxygen is part of a carbonyl (double bond to a carbon)
                # For alcohol fragment: oxygen must be terminal (exactly one carbon neighbor)
                nbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
                if role == 'alcohol':
                    if len(nbrs) != 1:
                        return False
                elif role == 'acid':
                    # Check if this oxygen is double-bonded to a carbon
                    is_carbonyl = False
                    for bond in atom.GetBonds():
                        # bond type check: bond.GetBondType() returns SINGLE, DOUBLE, etc.
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            is_carbonyl = True
                    if not is_carbonyl:
                        return False
        return True

    # Now fragment the molecule at the candidate ester bond.
    try:
        frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    except Exception as e:
        return False, f"Fragmentation failed at bond {bond_idx}: {str(e)}"
    
    fragments = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if len(fragments) != 2:
        return False, f"Fragmentation yielded {len(fragments)} fragments (expected 2)."
    
    # Identify fatty acid fragment (contains the carbonyl group) and fatty alcohol fragment.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)")
    
    frag_acid = None
    frag_alcohol = None
    for frag in fragments:
        if frag.HasSubstructMatch(acid_pattern):
            frag_acid = frag
        else:
            frag_alcohol = frag
    # If identification fails, try a fallback based on dummy atoms.
    if frag_acid is None or frag_alcohol is None:
        frag_acid, frag_alcohol = fragments
    
    # Check that fragments contain only C and O.
    if not is_fragment_aliphatic(frag_acid):
        return False, "Fatty acid fragment contains atoms other than C and O."
    if not is_fragment_aliphatic(frag_alcohol):
        return False, "Fatty alcohol fragment contains atoms other than C and O."
    
    # Both fragments must be acyclic.
    if not is_fragment_acyclic(frag_acid):
        return False, "Fatty acid fragment is cyclic; expected a linear aliphatic chain."
    if not is_fragment_acyclic(frag_alcohol):
        return False, "Fatty alcohol fragment is cyclic; expected a linear aliphatic chain."
    
    # Count aliphatic carbons.
    acid_carbon_count = count_aliphatic_carbons(frag_acid)
    alcohol_carbon_count = count_aliphatic_carbons(frag_alcohol)
    
    # Minimum required aliphatic carbons.
    min_acid_carbons = 12
    min_alcohol_carbons = 8
    if acid_carbon_count < min_acid_carbons:
        return False, f"Fatty acid fragment too short ({acid_carbon_count} aliphatic C; need >= {min_acid_carbons})."
    if alcohol_carbon_count < min_alcohol_carbons:
        return False, f"Fatty alcohol fragment too short ({alcohol_carbon_count} aliphatic C; need >= {min_alcohol_carbons})."
    
    # Enforce that each fragment is an unfunctionalized, linear chain.
    if not is_linear_aliphatic(frag_acid, "acid"):
        return False, "Fatty acid fragment is branched or functionally substituted beyond a normal acid chain."
    if not is_linear_aliphatic(frag_alcohol, "alcohol"):
        return False, "Fatty alcohol fragment is branched or functionally substituted beyond a normal alcohol chain."
    
    # Optionally, check minimal overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)."
    
    # If all tests pass, we declare the molecule as a wax ester.
    return True, (f"Candidate ester bond (bond idx {bond_idx}) splits molecule into a fatty acid fragment "
                  f"with {acid_carbon_count} aliphatic C and a fatty alcohol fragment with {alcohol_carbon_count} aliphatic C.")

# Example usage (uncomment to test):
# result, reason = is_wax_ester("CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC")
# print(result, reason)