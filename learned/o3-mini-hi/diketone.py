"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A diketone – a compound that contains exactly two independent ketone functionalities.
A “ketone functionality” here is defined as a carbonyl group R–CO–R where:
 • The carbonyl carbon is not aromatic (or part of an entirely aromatic ring);
 • It is sp2‐hybridized with exactly three neighbor atoms (one double‐bonded oxygen, two single-bonded carbons);
 • It has no attached hydrogens (to avoid picking up aldehydes).

In addition, if the two candidate ketone carbons lie together in a common ring
with mostly aromatic atoms (a quinone‐like conjugated system), we assume they are not the desired diketone functionalities.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone (contains exactly two independent ketone functionalities)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule has exactly two (non-quinone) ketone functionalities,
              False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring given as a tuple of atom indices).
    rings = mol.GetRingInfo().AtomRings()
    
    # Helper: returns True if an atom lies in a ring that is fully aromatic.
    def in_fully_aromatic_ring(atom):
        idx = atom.GetIdx()
        for ring in rings:
            if idx in ring:
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    return True
        return False

    # Define candidate ketone functionality:
    # Must be a carbon atom that is:
    #   - sp2-hybridized;
    #   - has exactly 3 neighbors (degree 3);
    #   - has no attached hydrogens;
    #   - is not aromatic and not part of a fully aromatic ring;
    #   - has exactly one double bond to an oxygen atom and two single bonds to carbon atoms.
    def is_candidate_ketone(atom):
        if atom.GetSymbol() != "C":
            return False
        # Check hybridization and degree.
        if atom.GetHybridization() != rdchem.HybridizationType.SP2:
            return False
        if atom.GetDegree() != 3:
            return False
        if atom.GetTotalNumHs() != 0:
            return False
        # Exclude if the atom itself is marked aromatic or is in a fully aromatic ring.
        if atom.GetIsAromatic() or in_fully_aromatic_ring(atom):
            return False

        doubleO_count = 0
        singleC_count = 0
        for bond in atom.GetBonds():
            btype = bond.GetBondType()
            nbr = bond.GetOtherAtom(atom)
            if btype == rdchem.BondType.DOUBLE:
                # Count only if neighbor is oxygen.
                if nbr.GetAtomicNum() == 8:
                    doubleO_count += 1
            elif btype == rdchem.BondType.SINGLE:
                if nbr.GetAtomicNum() == 6:
                    singleC_count += 1
        if doubleO_count == 1 and singleC_count == 2:
            return True
        return False

    # Gather candidate ketone atoms (store indices and atoms)
    candidate_atoms = []
    for atom in mol.GetAtoms():
        if is_candidate_ketone(atom):
            candidate_atoms.append(atom)
    
    ketone_count = len(candidate_atoms)
    
    # If not exactly 2 candidates, classification fails.
    if ketone_count != 2:
        return False, f"Found {ketone_count} ketone functionalities; exactly 2 are required for a diketone"
    
    # Extra check to avoid quinone-like systems:
    # If both candidate ketone atoms share at least one common ring where most atoms are aromatic,
    # then we assume they are part of a quinone-like conjugated system.
    candidate_idxs = {atom.GetIdx() for atom in candidate_atoms}
    for ring in rings:
        # Check if both candidate atoms are in this ring.
        if candidate_idxs.issubset(ring):
            # Count how many atoms in the ring are aromatic.
            aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            # If the ring is large enough and most atoms (say >50%) are aromatic, treat it as quinone-like.
            if len(ring) >= 5 and aromatic_count/len(ring) > 0.5:
                return False, "Ketone functionalities appear as part of a quinone-like aromatic system"
    
    return True, "Contains exactly two ketone functionalities"

# Example usage: run some tests if executed as a main program.
if __name__ == "__main__":
    test_cases = [
        ("O=C(CCCCCCCC)C(=O)C", "2,3-Undecanedione"),
        ("O=C(C(CCCCCCCC)CC)CC(=O)CCC", "7-Ethylpentadecane-4,6-dione"),
        ("CCCC(C(C)=O)=O", "2,3-Heptanedione"),
        ("O=C1C2CC(=O)C1(C)C2", "bornane-2,6-dione"),
        ("OC(=O)C(=O)CC(=O)C(O)=O", "2,4-dioxopentanedioic acid"),
    ]
    
    for s, name in test_cases:
        res, reason = is_diketone(s)
        print(f"SMILES: {s}\nNAME: {name}\nResult: {res}\nReason: {reason}\n{'-'*60}")