"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Any 1-benzopyran with an aryl substituent at position 3 (isoflavonoid).
An isoflavonoid is defined as a benzopyran (fused six-membered pyran ring with exactly one oxygen 
fused with a benzene ring sharing exactly 2–3 atoms) with an exocyclic aryl substituent attached 
to a non-fused pyran atom.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, otherwise False.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information (each ring is given as a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Identify candidate benzene rings:
    # These must be 6-membered, all atoms aromatic and all carbons.
    benzene_rings = []
    for ring in ring_info:
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() and 
                   mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
                   for idx in ring):
                benzene_rings.append(set(ring))
    
    if not benzene_rings:
        return False, "No benzene rings (aromatic 6-membered all-carbon) found for a fused system"
    
    # Identify candidate pyran rings:
    # These must be 6-membered rings with exactly one oxygen (atomic num 8) and the remainder carbons.
    pyran_rings = []
    for ring in ring_info:
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ocount = sum(1 for a in atoms if a.GetAtomicNum() == 8)
            ccount = sum(1 for a in atoms if a.GetAtomicNum() == 6)
            if ocount == 1 and ccount == 5:
                pyran_rings.append(set(ring))
    
    if not pyran_rings:
        return False, "No candidate pyran ring (6-membered with exactly one oxygen) found"

    # Look for a fused benzopyran system.
    # We require that one pyran and one benzene ring share 2 or 3 atoms (typical fusion in flavonoids).
    for pyran in pyran_rings:
        for benzene in benzene_rings:
            shared_atoms = pyran.intersection(benzene)
            if len(shared_atoms) in (2, 3):
                # We have a fused benzopyran system.
                fused_core = pyran.union(benzene)
                # The isoflavonoid motif requires an exocyclic aryl substituent attached
                # to the pyran ring at an atom that is outside the fused interface.
                non_fused_pyran_atoms = pyran.difference(shared_atoms)
                for atom_idx in non_fused_pyran_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    # Look at neighbors that are not in the pyran
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in pyran:
                            continue
                        # Neighbor must be aromatic and carbon-based to be a candidate aryl substituent.
                        if nbr.GetAtomicNum() != 6 or not nbr.GetIsAromatic():
                            continue
                        # Now, check that this neighbor is part of a 6-membered aromatic ring (i.e. resembles a phenyl)
                        nbr_found = False
                        for ring in ring_info:
                            if len(ring) == 6 and nbr.GetIdx() in ring:
                                if all(mol.GetAtomWithIdx(i).GetIsAromatic() and 
                                       mol.GetAtomWithIdx(i).GetAtomicNum() == 6 
                                       for i in ring):
                                    # Also check that this aromatic ring is not part of the fused core.
                                    if not set(ring).issubset(fused_core):
                                        nbr_found = True
                                        break
                        if nbr_found:
                            return True, ("Molecule contains a fused benzopyran core (a six-membered pyran with one oxygen fused with "
                                          "a benzene ring sharing 2–3 atoms) with an exocyclic aryl substituent.")
    # If no match is found, try to provide a more detailed reason.
    return False, "Scaffold not recognized as isoflavonoid (no fused benzopyran core with suitable exocyclic aryl substituent found)."

# Example usage (you may test with several SMILES strings):
if __name__ == '__main__':
    test_smiles = [
        # A couple of simple tests:
        "O1C(C2=CC=CC=C2)=CC(O)=C1",  # A simple benzopyran; may or may not have an extra aryl substituent.
        "COc1ccc(-c2coc3ccccc3c2)cc1" # Rough example: benzopyran plus a phenyl substituent.
    ]
    for smi in test_smiles:
        res, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")