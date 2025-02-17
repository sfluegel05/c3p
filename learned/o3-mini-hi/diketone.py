"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains exactly two ketone functionalities (diketone)
according to the definition: a “ketone” is a carbonyl group R–CO–R where the carbonyl 
carbon is not aromatic (or part of a fully aromatic ring), is trigonal (sp2, degree 3) with 
exactly one double bond to oxygen, no attached hydrogen (avoiding aldehydes), and is connected 
by two single bonds to carbon atoms. In addition, if a candidate ketone carbon is directly bonded 
to another candidate ketone carbon, we assume a quinone‐like conjugated system and ignore it.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone (contains exactly two independent ketone functionalities)
    based on its SMILES string.
    
    A ketone functionality (R–CO–R) is detected by first iterating over carbon atoms
    (ignoring those that are aromatic or part of a fully aromatic ring) and then requiring:
      - The atom is sp2 hybridized and has degree 3 (one double, two single connections)
      - It has no attached hydrogens (to avoid picking up aldehydes)
      - It is double bonded to exactly one oxygen atom (bond type DOUBLE)
      - Its two single-bonded neighbors are carbon atoms.
    In order to avoid flagging quinone-like systems (where two carbonyls are directly bonded
    to each other in a conjugated ring), we further check the single-bonded carbon neighbors.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule has exactly two ketone functionalities; False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ketone_count = 0
    
    # Get ring information once.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Helper: decide if an atom is "fully aromatic" in any ring. (If the atom lies in a ring that is 
    # entirely aromatic, then we wish to ignore the ketone.)
    def in_fully_aromatic_ring(atom):
        idx = atom.GetIdx()
        for ring in ring_info:
            if idx in ring:
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    return True
        return False

    # Helper: checks if a given atom qualifies as a ketone carbon by our criteria.
    def is_candidate_ketone(atom):
        # Must be carbon.
        if atom.GetSymbol() != "C":
            return False
        # Must be sp2 hybridized (trigonal) and degree 3.
        if atom.GetHybridization() != rdchem.HybridizationType.SP2:
            return False
        if atom.GetDegree() != 3:
            return False
        # Do not allow any attached hydrogens (avoid aldehydes, e.g.)
        if atom.GetTotalNumHs() != 0:
            return False
        # Skip if the atom is flagged as aromatic or in a fully aromatic ring.
        if atom.GetIsAromatic() or in_fully_aromatic_ring(atom):
            return False

        bonds = atom.GetBonds()
        doubleO = 0
        singleC = 0
        # Check bond details.
        for bond in bonds:
            btype = bond.GetBondType()
            nbr = bond.GetOtherAtom(atom)
            if btype == rdchem.BondType.DOUBLE:
                # Count only if neighbor is oxygen.
                if nbr.GetAtomicNum() == 8:
                    doubleO += 1
            elif btype == rdchem.BondType.SINGLE:
                if nbr.GetAtomicNum() == 6:
                    singleC += 1
        if doubleO == 1 and singleC == 2:
            return True
        return False

    # Iterate over atoms in the molecule.
    for atom in mol.GetAtoms():
        if is_candidate_ketone(atom):
            # Check for a quinone-like situation: if any of the two carbon neighbors (from single bonds)
            # is itself a candidate ketone carbon then we assume that this C=O is part of
            # a quinone-like conjugated system and ignore it.
            ignore_due_to_conjugation = False
            for bond in atom.GetBonds():
                if bond.GetBondType() == rdchem.BondType.SINGLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 6:
                        if is_candidate_ketone(nbr):
                            ignore_due_to_conjugation = True
                            break
            if ignore_due_to_conjugation:
                continue
            ketone_count += 1

    if ketone_count == 2:
        return True, "Contains exactly two ketone functionalities"
    else:
        return False, f"Found {ketone_count} ketone functionalities; exactly 2 are required for a diketone"

# Example usage: you can test a diketone SMILES.
if __name__ == "__main__":
    # For example, use 2,3-Undecanedione with SMILES "O=C(CCCCCCCC)C(=O)C"
    test_smiles = "O=C(CCCCCCCC)C(=O)C"
    result, reason = is_diketone(test_smiles)
    print(f"Result: {result}\nReason: {reason}")