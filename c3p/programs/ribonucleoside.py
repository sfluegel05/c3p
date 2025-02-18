"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: any nucleoside where the sugar component is D-ribose (ribonucleoside)

Definition:
A ribonucleoside is a nucleoside that has its nucleobase linked via an N-glycosidic bond to a D-ribofuranose unit.
The ribofuranose is a five-membered ring (4 carbons and 1 oxygen) and possesses an exocyclic CH2OH group 
attached to one of its carbons. Additionally, one ring carbon must be directly bonded to an aromatic nitrogen 
(as part of the nucleobase).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (a nucleoside where the sugar is D-ribose) 
    by checking for:
       - Absence of phosphorus (to avoid nucleotides).
       - Presence of a five-membered ring with exactly one oxygen and four carbons (ribofuranose).
       - At least one ring carbon having an exocyclic CH2OH group.
       - At least one ring carbon attached outside the ring to an aromatic nitrogen (nucleobase linkage).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ribonucleoside, else False.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Reject any molecules that contain phosphorus (atomic number 15)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Molecule contains phosphorus; likely a nucleotide rather than a nucleoside"

    # For better substituent perception, add explicit hydrogens.
    molH = Chem.AddHs(mol)
    ring_info = molH.GetRingInfo()

    # Define a helper function to decide if an exocyclic atom is a CH2OH substituent.
    def is_exocyclic_CH2OH(atom):
        # Must be a carbon atom.
        if atom.GetAtomicNum() != 6:
            return False
        # Check that it has exactly two hydrogens.
        if atom.GetTotalNumHs() != 2:
            return False
        # It must be attached to an oxygen that bears at least one hydrogen (i.e. is -OH).
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                return True
        return False

    # Iterate over each ring in the molecule to find candidate ribofuranose rings.
    for ring in ring_info.AtomRings():
        # Consider only rings of exactly five atoms.
        if len(ring) != 5:
            continue

        # Get the atoms of the ring.
        ring_atoms = [molH.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen and carbon atoms within the ring.
        num_oxygens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        # The ring must contain exactly 1 oxygen and 4 carbons.
        if num_oxygens != 1 or num_carbons != 4:
            continue

        # Two features we need for ribose:
        #  - exocyclic CH2OH group attached to one ring atom
        #  - a linkage to a nucleobase as indicated by an aromatic nitrogen attached to a ring atom
        ch2oh_found = False
        base_link_found = False

        # For each atom in this candidate ring...
        for idx in ring:
            ring_atom = molH.GetAtomWithIdx(idx)
            # Get neighbors outside of the ring.
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Skip atoms that are part of the ring.
                # Check for CH2OH substituent.
                if is_exocyclic_CH2OH(nbr):
                    ch2oh_found = True
                # Check for nucleobase linkage: an aromatic nitrogen.
                if nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                    base_link_found = True

        # If both substituents are detected for this ring, we have found a ribofuranose unit.
        if ch2oh_found and base_link_found:
            return True, ("Found a ribofuranose ring (1 O and 4 C) with an exocyclic CH2OH group and "
                          "a nucleobase linkage via an aromatic nitrogen")
    
    return False, "No ribose moiety (with CH2OH and nucleobase linkage) found"


# Example usage and testing (if run as __main__)
if __name__ == "__main__":
    examples = {
        "1-methyladenosine": "Cn1cnc2n(cnc2c1=N)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
        "3,4-dihydrozebularine": "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=CCNC1=O",
        "nucleocidin": "NC1=C2N=CN([C@@H]3O[C@](F)(COS(N)(=O)=O)[C@@H](O)[C@H]3O)C2=NC=N1",
        "cytidine": "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1",
        "aminodeoxyfutalosine": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CCC(=O)c2cccc(c2)C(O)=O)[C@@H](O)[C@H]1O",
        "5'-deoxyadenosine (should fail ribonucleoside)": "C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
        "Adenosylcobinamide phosphate (should be rejected)": "C[C@H](CNC(=O)CC[C@]1(C)[C@@H](CC(N)=O)[C@H]2N3C1=C(C)C1=[N+]4C(=CC5=[N+]6C(=C(C)C7=[N+]([C@]2(C)[C@@](C)(CC(N)=O)[C@@H]7CCC(N)=O)[Co--]346C[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@](C)(CC(N)=O)[C@@H]5CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O)OP(O)(O)=O"
    }
    
    for name, smi in examples.items():
        is_ribo, reason = is_ribonucleoside(smi)
        print(f"{name}: {is_ribo} ({reason})")