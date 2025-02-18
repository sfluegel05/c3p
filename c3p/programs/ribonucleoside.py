"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: any nucleoside where the sugar component is D-ribose (ribonucleoside)
Definition: A ribonucleoside is a nucleoside having its nucleobase linked
             via an N-glycosidic bond to a D-ribofuranose (a five-membered ring
             containing 1 oxygen and 4 carbons, plus an exocyclic CH2OH group at the 5'-position).
Improvements:
  1. Reject compounds containing phosphorus.
  2. Look for five-membered rings that have exactly 1 oxygen and 4 carbons.
  3. For each such ring, check that at least one ring carbon has an exocyclic CH2OH substituent.
  4. Also require that at least one ring carbon is directly bonded (outside the ring) to
     an aromatic nitrogen – taken as a proxy for a nucleobase linkage.
If all these conditions are met, then the molecule is classified as a ribonucleoside.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (nucleoside with D-ribose sugar)
    by checking for (i) the absence of phosphorus, (ii) the presence of a five-membered ring
    with exactly one oxygen and four carbons, (iii) existence of an exocyclic CH2OH group (the 5'-substituent),
    and (iv) attachment of the sugar ring to an aromatic nitrogen (i.e. nucleobase linkage).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a ribonucleoside, else False.
        str: Explanation for the classification decision.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Reject molecules that contain phosphorus (to avoid nucleotides)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely a nucleotide rather than a nucleoside"
    
    # For better perception of substituents, add hydrogens.
    molH = Chem.AddHs(mol)
    ring_info = molH.GetRingInfo()
    
    # Helper: Check if a candidate exocyclic substituent is a CH2OH group.
    def is_CH2OH(atom):
        # We expect atom to be carbon (sp3) with exactly two hydrogens [CH2] and an -OH group.
        if atom.GetAtomicNum() != 6:
            return False
        # Count the total number of hydrogens explicitly present.
        # (Using GetTotalNumHs() gives the total implicit+explicit H count.)
        if atom.GetTotalNumHs() < 2:
            return False
        # Now, look for an oxygen neighbor (outside the sugar ring) that itself has at least one hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check if this oxygen is an OH (attached to at least one hydrogen)
                if nbr.GetTotalNumHs() >= 1:
                    return True
        return False
    
    # For each five-membered ring, try to detect a ribofuranose unit.
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        # Count oxygen and carbon atoms in this ring.
        ring_atoms = [molH.GetAtomWithIdx(idx) for idx in ring]
        oxygens = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        carbons = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]
        if len(oxygens) != 1 or len(carbons) != 4:
            continue  # not a furanose pattern
        
        # Check for the typical 5'-CH2OH substituent.
        # In a D-ribofuranose, one of the ring carbons (usually C4') is attached outside the ring
        # to a CH2OH group.
        ch2oh_found = False
        for idx in ring:
            atom = molH.GetAtomWithIdx(idx)
            # Look at neighbors that are NOT part of the ring:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # We want a carbon that qualifies as CH2OH.
                if is_CH2OH(nbr):
                    ch2oh_found = True
                    break
            if ch2oh_found:
                break
        if not ch2oh_found:
            # This ring does not display the CH2OH substituent; skip it.
            continue
        
        # Check for attachment of the sugar ring to a nucleobase:
        # Look for any ring carbon (candidate for the anomeric carbon) that
        # is bonded (outside the ring) to an aromatic nitrogen.
        nucleobase_attached = False
        for idx in ring:
            atom = molH.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only consider carbons
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                    nucleobase_attached = True
                    break
            if nucleobase_attached:
                break
        
        if nucleobase_attached:
            return True, ("Found five-membered ribofuranose ring (1 O, 4 C) with a CH2OH substituent "
                          "and a ring carbon attached to an aromatic nitrogen (nucleobase linkage)")
    
    return False, "No ribose moiety (with CH2OH and nucleobase linkage) found"

# Example usage and testing:
if __name__ == "__main__":
    # Example SMILES strings and expected outcomes (per provided list)
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
        res, reason = is_ribonucleoside(smi)
        print(f"{name}: {res} ({reason})")