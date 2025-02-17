"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: any nucleoside where the sugar component is D-ribose (ribonucleoside)

Definition:
A ribonucleoside is a nucleoside that has its nucleobase linked via an N-glycosidic
bond to a D-ribofuranose unit. The ribofuranose is a five-membered ring (4 carbons and 1 oxygen)
and possesses an exocyclic CH2OH group attached to one of the ring carbons (typically the 5'-position).
Additionally, one ring carbon must be directly bonded to an aromatic nitrogen (as part of the nucleobase).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (a nucleoside where the sugar is D-ribose)
    by checking for:
       - Absence of phosphorus (to exclude nucleotides).
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

    # Add explicit hydrogens to help perceiving substituents.
    molH = Chem.AddHs(mol)
    ring_info = molH.GetRingInfo()

    # Helper function to check if an atom (a candidate substituent) is a CH2OH group.
    def is_CH2OH_group(atom):
        # The candidate must be a carbon.
        if atom.GetAtomicNum() != 6:
            return False
        # For a CH2 group, expect exactly 2 hydrogens.
        if atom.GetTotalNumHs() != 2:
            return False
        # Now, one of its neighbors (besides the connecting ring atom) should be an oxygen
        # that carries at least one hydrogen (i.e. an -OH group).
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                if nbr.GetTotalNumHs() >= 1:
                    return True
        return False

    # Now, search for a 5-membered ring that fits the criteria for ribofuranose.
    for ring in ring_info.AtomRings():
        # consider only rings of exactly five atoms
        if len(ring) != 5:
            continue
        # Get the atoms in the ring.
        ring_atoms = [molH.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen and carbon atoms in the ring.
        num_oxygens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        # For a ribofuranose, we need exactly 1 oxygen and 4 carbons.
        if num_oxygens != 1 or num_carbons != 4:
            continue

        # Features needed for ribose: an exocyclic CH2OH and a nucleobase linkage.
        ch2oh_found = False
        base_link_found = False
        
        # Iterate each atom in the candidate ring.
        for idx in ring:
            ring_atom = molH.GetAtomWithIdx(idx)
            # Check all neighbors outside the ring.
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # ignore atoms that are in the ring
                # Check for exocyclic CH2OH substituent.
                if not ch2oh_found and is_CH2OH_group(nbr):
                    ch2oh_found = True
                # Check for nucleobase linkage: an aromatic nitrogen.
                # We assume that a nucleobase will be attached via an N-glycosidic bond.
                if not base_link_found and nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                    base_link_found = True
                    
        if ch2oh_found and base_link_found:
            return True, ("Found a ribofuranose ring (1 O and 4 C) with an exocyclic CH2OH group and "
                          "a nucleobase linkage via aromatic nitrogen")
    
    return False, "No ribose moiety (with CH2OH and nucleobase linkage) found"

    
# Example usage and testing:
if __name__ == "__main__":
    examples = {
        "1-methyladenosine": "Cn1cnc2n(cnc2c1=N)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
        "3,4-dihydrozebularine": "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=CCNC1=O",
        "nucleocidin": "NC1=C2N=CN([C@@H]3O[C@](F)(COS(N)(=O)=O)[C@@H](O)[C@H]3O)C2=NC=N1",
        "cytidine": "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1",
        "aminodeoxyfutalosine": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CCC(=O)c2cccc(c2)C(O)=O)[C@@H](O)[C@H]1O",
        "5'-deoxyadenosine (should be rejected)": "C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"
    }
    
    for name, smi in examples.items():
        is_ribo, reason = is_ribonucleoside(smi)
        print(f"{name}: {is_ribo} ({reason})")