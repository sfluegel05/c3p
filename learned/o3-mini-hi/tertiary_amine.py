"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: Tertiary amine
A tertiary amine is defined as a compound derived formally from ammonia by replacing all three hydrogens 
with hydrocarbyl groups. In our code the molecule must have at least one neutral (formal charge 0) 
sp³ nitrogen atom that has exactly three heavy (non‐H) neighbors – all carbons – and no attached hydrogens.
Moreover, if the candidate nitrogen is not in a ring then removing it should disconnect its three substituents.
This extra “isolation” test reduces the risk that a structure containing a tertiary amine sub‐motif 
within a larger scaffold is (mis)classified as “tertiary amine”.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines if the molecule (given as a SMILES string) is a tertiary amine type compound.
    A compound is considered a tertiary amine if it is formally derived from ammonia
    by replacing all three hydrogens with hydrocarbyl groups. In our implementation we
    look for a neutral sp3 nitrogen that has exactly three heavy (non‐hydrogen) neighbors – all carbons –
    none of which appears to be a carbonyl carbon (C double‐bonded to O). In addition, if the candidate is not
    in a ring then we require that removing the nitrogen disconnects the substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a tertiary amine (by our criteria), False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES and add explicit hydrogens to get correct counts.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol_with_H = Chem.AddHs(mol)
    
    # Candidate tertiary amine atom will be stored here.
    candidate = None

    # Iterate over atoms looking for a candidate nitrogen.
    for atom in mol_with_H.GetAtoms():
        if atom.GetAtomicNum() != 7:  # must be nitrogen
            continue
        if atom.GetFormalCharge() != 0:
            continue
        # Only consider sp3 atoms (this rule will reject most aromatic or planar N’s)
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # Check that there is no hydrogen attached (explicitly or implicitly)
        # (Since we added Hs, GetDegree() counts all bonds)
        # We require that none of the neighbors is hydrogen.
        if any(neigh.GetAtomicNum() == 1 for neigh in atom.GetNeighbors()):
            continue
        # Now, count heavy (non-H) neighbors.
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 3:
            continue
        # Require that all three neighbors are carbon.
        if not all(n.GetAtomicNum() == 6 for n in heavy_neighbors):
            continue
        # Exclude cases where one of the carbon neighbors is a carbonyl carbon
        # i.e. it has a double bond to an oxygen.
        is_carbonyl = False
        for nbr in heavy_neighbors:
            for bond in nbr.GetBonds():
                # Check for a double bond
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8:
                        is_carbonyl = True
                        break
            if is_carbonyl:
                break
        if is_carbonyl:
            continue
        
        # We have a candidate tertiary amine nitrogen.
        candidate = atom
        break

    if candidate is None:
        return False, "No tertiary amine group found in the molecule."
    
    # For molecules where the candidate N is not in a ring,
    # we require that its three substituents are only connected via the candidate.
    if not candidate.IsInRing():
        # Make a copy in which we remove the candidate nitrogen.
        emol = Chem.EditableMol(mol_with_H)
        cand_idx = candidate.GetIdx()
        emol.RemoveAtom(cand_idx)
        mol_removed = emol.GetMol()
        # Get connected fragments (each returned as a tuple of atom indices)
        frags = Chem.GetMolFrags(mol_removed, asMols=False)
        # Get the indices (in mol_with_H) of the candidate's original neighbors.
        nbr_indices = sorted(n.GetIdx() for n in candidate.GetNeighbors())
        # In the molecule without the candidate, each of these three neighbors should lie in a different fragment.
        frag_assignment = {}  # atom index -> fragment id
        for frag_id, frag in enumerate(frags):
            for idx in frag:
                frag_assignment[idx] = frag_id
        # For each neighbor, get its fragment.
        nbr_frags = [frag_assignment.get(idx, -1) for idx in nbr_indices]
        if len(set(nbr_frags)) != 3:
            return False, ("Tertiary amine found at atom index {} appears integrated into a larger "
                           "scaffold (upon removal its substituents are not fully disconnected)."
                           .format(candidate.GetIdx()))
    # If we got here the candidate N qualifies.
    return True, ("Found tertiary amine at atom index {} with three carbon substituents."
                  .format(candidate.GetIdx()))

# Example usage:
if __name__ == "__main__":
    # Try a few examples:
    examples = [
        ("Tri-allate", "CC(C)N(C(C)C)C(=O)SCC(Cl)=C(Cl)Cl"),
        ("triethylamine", "CCN(CC)CC"),
        ("(R)-fenpropidin", "C[C@@H](CN1CCCCC1)Cc1ccc(cc1)C(C)(C)C"),
        ("3-quinuclidinol", "OC1C[N@@]2CC[C@H]1CC2"),
        ("N,N-dimethylethanolamine", "CN(C)CCO"),
        # A false negative test shows no tertiary amine:
        ("Azoxymycin A", "O=C(NC(C(=O)O)CCC(=O)N)/C=C/C=C/C1=CC=C([N+]([O-])=NC2=CC=C(/C=C/C=C/C(=O)NC(C(=O)O)CCC(=O)N)C=C2)C=C1")
    ]
    for name, smi in examples:
        result, reason = is_tertiary_amine(smi)
        print(f"{name}: {result} => {reason}")